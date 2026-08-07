#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(HDF5Array)
  library(GenomicRanges)
  library(yaml)
})

ROOT <- "/projects/b1213/users/kynon/projects/dna-methylation-heritability"
source(file.path(ROOT, "inputs/cell_proportions/_h/dnam_deconvolution_utils.R"))
args <- commandArgs(trailingOnly = TRUE)
reg <- normalize_region(if (length(args)) args[[1]] else "caudate")
signature_platform <- if (length(args) >= 2L) toupper(args[[2]]) else "850K"
use_coordinate_wgbs <- signature_platform == "WGBS"
cfg <- yaml::read_yaml(file.path(ROOT, "config/cell_deconvolution.yml"))
if (!reg %in% cfg$regions) stop("Unsupported region: ", reg)

render <- function(template, region) gsub("\\{region\\}", region, template)
ref_dir <- file.path(ROOT, cfg$paths$reference_dir)
work_template <- if (use_coordinate_wgbs) cfg$paths$dnam_wgbs_work_template else cfg$paths$dnam_work_template
work_dir <- file.path(ROOT, render(work_template, reg))
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

marker_file <- file.path(ref_dir, if (use_coordinate_wgbs) {
  "scmd_coordinate_wgbs_marker_map.tsv"
} else {
  "scmd_wgbs_marker_map.tsv"
})
if (!file.exists(marker_file)) stop("Missing marker map; run 04.prepare_scmd_markers.R")
marker_map <- fread(marker_file)
marker_map <- marker_map[eligible == TRUE]
if (!nrow(marker_map)) stop("No eligible reference markers")

phen <- fread(file.path(ROOT, cfg$paths$phenotype_table))
phen[, region_norm := normalize_region(region)]
eligible_samples <- unique(as.character(phen[region_norm == reg]$brnum))

h5_path <- file.path(ROOT, render(cfg$paths$hdf5_bsseq_template, reg))
if (!dir.exists(h5_path)) stop("Missing HDF5-backed BSseq object: ", h5_path)
bs <- loadHDF5SummarizedExperiment(h5_path)
if (!"brnum" %in% colnames(colData(bs))) stop("BSseq colData lacks brnum")
bs_ids <- as.character(colData(bs)$brnum)
sample_ids <- eligible_samples[eligible_samples %in% bs_ids]
if (!length(sample_ids)) stop("No phenotype/BSseq sample overlap for ", reg)
bs_col <- match(sample_ids, bs_ids)

coord <- unique(marker_map[, .(seqnames_hg38, pos_hg38)])
coord[, coord_key := paste(seqnames_hg38, pos_hg38, sep = ":")]
query <- GRanges(coord$seqnames_hg38, IRanges(coord$pos_hg38, width = 1))
hit <- findOverlaps(query, rowRanges(bs), type = "equal", select = "first", ignore.strand = TRUE)
coord[, bs_row := as.integer(hit)]
coord[, found_in_bsseq := !is.na(bs_row)]

found <- coord[found_in_bsseq == TRUE]
if (!nrow(found)) stop("No lifted scMD markers overlap ", reg, " BSseq rows")
m <- as.matrix(assay(bs, "M")[found$bs_row, bs_col, drop = FALSE])
cov <- as.matrix(assay(bs, "Cov")[found$bs_row, bs_col, drop = FALSE])
colnames(m) <- colnames(cov) <- sample_ids
rownames(m) <- rownames(cov) <- found$coord_key

min_cov <- as.numeric(cfg$filtering$min_coverage)
usable <- is.finite(cov) & cov >= min_cov
marker_fraction <- rowMeans(usable)
marker_pass <- marker_fraction >= as.numeric(cfg$filtering$min_marker_sample_fraction)
found[, `:=`(
  fraction_samples_cov_ge_min = marker_fraction,
  mean_coverage = rowMeans(cov, na.rm = TRUE),
  marker_coverage_pass = marker_pass
)]

pass_rows <- which(marker_pass)
sample_fraction <- colMeans(usable[pass_rows, , drop = FALSE])
sample_pass <- sample_fraction >= as.numeric(cfg$filtering$min_sample_marker_fraction)
sample_qc <- data.table(
  sample_id = sample_ids,
  region = reg,
  n_marker_coordinates = length(pass_rows),
  n_marker_coordinates_covered = colSums(usable[pass_rows, , drop = FALSE]),
  marker_coverage_fraction = sample_fraction,
  mean_marker_coverage = colMeans(cov[pass_rows, , drop = FALSE], na.rm = TRUE),
  sample_qc_pass = sample_pass
)

keep_samples <- which(sample_pass)
if (!length(keep_samples)) stop("All samples failed marker coverage QC")
beta_coord <- m[pass_rows, keep_samples, drop = FALSE] / cov[pass_rows, keep_samples, drop = FALSE]
beta_coord[!usable[pass_rows, keep_samples, drop = FALSE]] <- NA_real_
for (i in seq_len(nrow(beta_coord))) {
  missing <- !is.finite(beta_coord[i, ])
  if (any(missing)) beta_coord[i, missing] <- median(beta_coord[i, !missing], na.rm = TRUE)
}

passed_coord <- found[pass_rows]
passed_coord[, coord_key := paste(seqnames_hg38, pos_hg38, sep = ":")]
target <- unique(marker_map[, .(reference, target_id, seqnames_hg38, pos_hg38)])
target[, coord_key := paste(seqnames_hg38, pos_hg38, sep = ":")]
target <- target[coord_key %in% passed_coord$coord_key]
target <- target[!duplicated(target_id)]
idx <- match(target$coord_key, rownames(beta_coord))
bulk <- beta_coord[idx, , drop = FALSE]
rownames(bulk) <- target$target_id
bulk <- bulk[!duplicated(rownames(bulk)), , drop = FALSE]

marker_qc <- merge(
  coord,
  found[, .(coord_key, fraction_samples_cov_ge_min, mean_coverage, marker_coverage_pass)],
  by = "coord_key", all.x = TRUE
)
fwrite(marker_qc, file.path(work_dir, "marker_coverage_qc.tsv"), sep = "\t", na = "NA")
fwrite(sample_qc, file.path(work_dir, "sample_marker_qc.tsv"), sep = "\t")
saveRDS(bulk, file.path(work_dir, "scmd_bulk_markers.rds"), compress = FALSE)
fwrite(as.data.table(bulk, keep.rownames = "target_id"),
       file.path(work_dir, "scmd_bulk_markers.tsv.gz"), sep = "\t")

summary <- data.table(
  region = reg,
  n_expected_samples = length(eligible_samples),
  n_bsseq_overlap_samples = length(sample_ids),
  n_samples_passing = sum(sample_pass),
  n_reference_candidate_coordinates = nrow(coord),
  n_coordinates_found = nrow(found),
  n_coordinates_passing_coverage = sum(marker_pass),
  n_target_ids_in_bulk = nrow(bulk),
  signature_platform = signature_platform,
  min_coverage = min_cov
)
fwrite(summary, file.path(work_dir, "extraction_summary.tsv"), sep = "\t")
capture.output(sessionInfo(), file = file.path(work_dir, "extraction_session_info.txt"))
print(summary)
