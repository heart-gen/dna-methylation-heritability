#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(minfi)
  library(rtracklayer)
  library(yaml)
})

ROOT <- "/projects/b1213/users/kynon/projects/dna-methylation-heritability"
source(file.path(ROOT, "inputs/cell_proportions/_h/dnam_deconvolution_utils.R"))
cfg <- yaml::read_yaml(file.path(ROOT, "config/cell_deconvolution.yml"))
ref_dir <- file.path(ROOT, cfg$paths$reference_dir)
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(file.path(ref_dir, "R_libs"), .libPaths()))
args <- commandArgs(trailingOnly = TRUE)
use_coordinate_wgbs <- length(args) && toupper(args[[1]]) == "WGBS"

required <- if (use_coordinate_wgbs) {
  c("Lee_7ct_WGBS.rda", "Tian_7ct_WGBS.rda")
} else {
  c("Lee_7ct_450850.rda", "Tian_7ct_450850.rda")
}
missing <- required[!file.exists(file.path(ref_dir, required))]
if (length(missing)) stop("Missing scMD reference assets: ", paste(missing, collapse = ", "))

refs <- if (use_coordinate_wgbs) {
  list(
    Lee = load_reference_objects(file.path(ref_dir, "Lee_7ct_WGBS.rda"), "Lee", "WGBS"),
    Tian = load_reference_objects(file.path(ref_dir, "Tian_7ct_WGBS.rda"), "Tian", "WGBS")
  )
} else {
  list(
    Lee = load_reference_objects(file.path(ref_dir, "Lee_7ct_450850.rda"), "Lee", "450850"),
    Tian = load_reference_objects(file.path(ref_dir, "Tian_7ct_450850.rda"), "Tian", "450850")
  )
}
candidate_n <- as.integer(cfg$reference$candidate_markers_per_cell)
markers <- rbindlist(lapply(names(refs), function(nm) {
  top_reference_markers(refs[[nm]], nm, candidate_n)
}))
if (use_coordinate_wgbs) {
  coords <- parse_reference_coordinates(markers$target_id)
  markers <- cbind(markers, coords)
  annotation_package <- NA_character_
  markers[, signature_platform := "WGBS"]
} else {
  annotation_package <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
  if (!requireNamespace(annotation_package, quietly = TRUE)) {
    annotation_package <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
    if (!requireNamespace(annotation_package, quietly = TRUE)) {
      stop("Neither the EPIC nor 450K hg19 minfi annotation package is installed")
    }
    warning("EPIC annotation unavailable; using the prespecified 450K coordinate fallback")
  }
  suppressPackageStartupMessages(library(annotation_package, character.only = TRUE))
  annotation_object <- get(annotation_package)
  anno <- minfi::getAnnotation(annotation_object)
  anno <- as.data.table(as.data.frame(anno), keep.rownames = "target_id")
  coords <- anno[, .(
    target_id,
    seqnames_hg19 = as.character(chr),
    pos_hg19 = as.integer(pos)
  )]
  markers <- merge(markers, coords, by = "target_id", all.x = TRUE, sort = FALSE)
  markers[, signature_platform := if (grepl("450kanno", annotation_package)) {
    as.character(cfg$reference$fallback_signature_platform)
  } else {
    as.character(cfg$reference$primary_signature_platform)
  }]
}

keys <- unique(markers[, .(reference, target_id, seqnames_hg19, pos_hg19)])
valid_coord <- !is.na(keys$pos_hg19) & grepl("^chr([0-9]+|X)$", keys$seqnames_hg19)
keys[, mapping_status := fifelse(valid_coord, "pending", "invalid_source_coordinate")]

materialize_chain <- function(path) {
  if (!grepl("\\.gz$", path)) return(path)
  out <- tempfile(fileext = ".chain")
  input <- gzfile(path, "rb")
  output <- file(out, "wb")
  on.exit({close(input); close(output)}, add = TRUE)
  repeat {
    chunk <- readBin(input, what = "raw", n = 1024L * 1024L)
    if (!length(chunk)) break
    writeBin(chunk, output)
  }
  out
}

eligible_source <- which(valid_coord)
source_gr <- GRanges(
  seqnames = keys$seqnames_hg19[eligible_source],
  ranges = IRanges(keys$pos_hg19[eligible_source], width = 1),
  key_index = eligible_source
)
chain_path <- materialize_chain(file.path(ROOT, cfg$paths$hg19_to_hg38_chain))
chain <- import.chain(chain_path)
lifted <- liftOver(source_gr, chain)
lift_n <- lengths(lifted)

keys[, `:=`(seqnames_hg38 = NA_character_, pos_hg38 = NA_integer_)]
unique_lift <- which(lift_n == 1L)
lifted_one <- unlist(lifted[unique_lift], use.names = FALSE)
dest_rows <- mcols(source_gr)$key_index[unique_lift]
keys[dest_rows, `:=`(
  seqnames_hg38 = as.character(seqnames(lifted_one)),
  pos_hg38 = as.integer(start(lifted_one)),
  mapping_status = "unique_liftover"
)]
keys[mcols(source_gr)$key_index[lift_n == 0L], mapping_status := "unmapped"]
keys[mcols(source_gr)$key_index[lift_n > 1L], mapping_status := "ambiguous_liftover"]

mapped <- which(keys$mapping_status == "unique_liftover")
dest_gr <- GRanges(
  keys$seqnames_hg38[mapped],
  IRanges(keys$pos_hg38[mapped], width = 1),
  key_index = mapped
)

overlap_flag <- function(query, path) {
  if (!file.exists(path)) stop("Missing exclusion asset: ", path)
  subject <- rtracklayer::import(path)
  overlapsAny(query, subject, ignore.strand = TRUE)
}

keys[, `:=`(overlaps_blacklist = FALSE, overlaps_ct_snp = FALSE)]
if (isTRUE(cfg$filtering$exclude_blacklist) && length(dest_gr)) {
  keys[mcols(dest_gr)$key_index, overlaps_blacklist := overlap_flag(
    dest_gr, file.path(ROOT, cfg$paths$blacklist_hg38)
  )]
}
if (isTRUE(cfg$filtering$exclude_ct_snps) && length(dest_gr)) {
  keys[mcols(dest_gr)$key_index, overlaps_ct_snp := overlap_flag(
    dest_gr, file.path(ROOT, cfg$paths$ct_snps_hg38)
  )]
}

keys[, duplicated_reference_coordinate := FALSE]
dup <- keys[mapping_status == "unique_liftover",
            .N, by = .(reference, seqnames_hg38, pos_hg38)][N > 1L]
if (nrow(dup)) {
  keys[dup, on = .(reference, seqnames_hg38, pos_hg38),
       duplicated_reference_coordinate := TRUE]
}
keys[, eligible := mapping_status == "unique_liftover" &
       !overlaps_blacklist & !overlaps_ct_snp & !duplicated_reference_coordinate]

markers <- merge(markers, keys, by = c("reference", "target_id", "seqnames_hg19", "pos_hg19"),
                 all.x = TRUE, sort = FALSE)
setorder(markers, reference, cell_type, marker_rank)

availability <- markers[, .(
  n_candidates = uniqueN(target_id),
  n_unique_liftover = uniqueN(target_id[mapping_status == "unique_liftover"]),
  n_blacklist = uniqueN(target_id[overlaps_blacklist]),
  n_ct_snp = uniqueN(target_id[overlaps_ct_snp]),
  n_duplicate_coordinate = uniqueN(target_id[duplicated_reference_coordinate]),
  n_eligible = uniqueN(target_id[eligible])
), by = .(reference, cell_type)]

file_tag <- if (use_coordinate_wgbs) "scmd_coordinate_wgbs" else "scmd_wgbs"
fwrite(markers, file.path(ref_dir, paste0(file_tag, "_marker_map.tsv")), sep = "\t", na = "NA")
fwrite(availability, file.path(ref_dir, paste0(file_tag, "_marker_availability.tsv")), sep = "\t")
fwrite(data.table(
  signature_platform = unique(markers$signature_platform),
  annotation_package = annotation_package,
  annotation_version = if (is.na(annotation_package)) NA_character_ else as.character(utils::packageVersion(annotation_package))
), file.path(ref_dir, if (use_coordinate_wgbs) {
  "scmd_coordinate_wgbs_signature_platform.tsv"
} else {
  "scmd_signature_platform.tsv"
}), sep = "\t")
capture.output(sessionInfo(), file = file.path(ref_dir, "marker_preparation_session_info.txt"))

if (any(availability$n_eligible < cfg$reference$markers_per_cell)) {
  warning("Some reference/cell-type pools have fewer than the requested final marker count; see availability table")
}
cat("Prepared", uniqueN(markers$target_id), "reference marker IDs;",
    uniqueN(markers[eligible == TRUE]$target_id), "passed coordinate/exclusion QC\n")
