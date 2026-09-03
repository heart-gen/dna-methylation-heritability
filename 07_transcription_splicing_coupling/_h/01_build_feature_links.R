#!/usr/bin/env Rscript
#### 07_transcription_splicing_coupling -- build the tested universe ####
##
## Usage:
##   Rscript _h/01_build_feature_links.R --run-id tsc-AA-caudate-YYYYMMDD
##
## Builds VMR-to-feature links against the ACCEPTED Module 01 VMR boundaries.
##
## The legacy analysis read a precomputed link table
## (`{tissue}_vmr_genes_within_250kb_hg38.tsv`) keyed to pre-repair VMRs. Those
## links cannot be reused: AGENTS.md 6 states that correcting VMR definitions
## "requires recomputing every analysis based on VMR membership, boundaries,
## summaries, or classification", and a link is a function of the boundary.
##
## Writes the tested universe explicitly (n VMRs, n features, n pairs, per
## modality), because AGENTS.md 7.6 requires the denominator of a coupling claim
## to be documented rather than implied.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE", file.path(Sys.getenv("V2_REPO_ROOT", "."), "07_transcription_splicing_coupling", "_h")), "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
})

MODULE <- "07_transcription_splicing_coupling"

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
## Select outside the data.table `[`: inside it, a bare `field` binds to the
## column, not the argument (see 00_shared/gates.R for the same trap).
mf <- function(field, required = TRUE) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0 || is.na(v[1])) {
        if (required) stop("Manifest has no value for '", field, "'")
        return(NA_character_)
    }
    v[1]
}

cohort <- mf("cohort"); region <- mf("region")
vmr_run <- mf("upstream_vmr_catalog_run_id")
ts <- load_run_config("transcription_splicing", run_dir)
enabled <- strsplit(mf("modalities"), ",", fixed = TRUE)[[1]]

norm_chr <- function(x) {
    x <- as.character(x)
    ifelse(startsWith(x, "chr"), x, paste0("chr", x))
}

## ------------------------------------------------------------------- VMRs
vmr_bed <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs", vmr_run,
                     "vmr", "vmr.bed")
if (!file.exists(vmr_bed)) stop("Accepted VMR bed not found: ", vmr_bed)
vmr <- fread(vmr_bed, header = FALSE, col.names = c("chrom", "start", "end"))
vmr[, chrom := norm_chr(chrom)]
## Canonical v2 VMR identifier, matching 01_vmr_catalog/_h/02_summarize.R and
## the IDs carried by Modules 02 and 05. The legacy underscore form
## (chr1_134100_134201) survives only as the .phen FILENAME convention, and
## 02_run_local_associations.R converts to it when reading those files.
vmr[, vmr_id := paste0(chrom, ":", start, "-", end)]
vmr <- unique(vmr, by = "vmr_id")
message("[07] ", nrow(vmr), " VMRs from ", vmr_run)

vmr_gr <- GRanges(vmr$chrom, IRanges(vmr$start + 1L, vmr$end))
mcols(vmr_gr)$vmr_id <- vmr$vmr_id

#' Local links between VMRs and annotated features within `window` bp.
#'
#' Distance is measured between intervals, so a feature overlapping the VMR has
#' distance 0. Both the link and its distance are kept: AGENTS.md 7.6 requires
#' VMR-to-feature distance as a covariate in the coupling model.
window_links <- function(feat, window) {
    feat <- feat[is.finite(start) & is.finite(end) & nzchar(chrom)]
    feat[, chrom := norm_chr(chrom)]
    fgr <- GRanges(feat$chrom, IRanges(pmin(feat$start, feat$end),
                                       pmax(feat$start, feat$end)))
    mcols(fgr)$feature_id <- feat$feature_id
    hits <- findOverlaps(vmr_gr, fgr, maxgap = as.integer(window))
    if (length(hits) == 0) return(data.table())
    dt <- data.table(
        vmr_id = mcols(vmr_gr)$vmr_id[queryHits(hits)],
        feature_id = mcols(fgr)$feature_id[subjectHits(hits)],
        distance = GenomicRanges::distance(vmr_gr[queryHits(hits)],
                                           fgr[subjectHits(hits)])
    )
    dt[is.na(distance), distance := 0L]
    unique(dt, by = c("vmr_id", "feature_id"))
}

gene_annot <- fread(file.path(repo_root(), ts$annotation$gene))
gene_annot[, feature_id := gene_id]
psi_annot <- fread(file.path(repo_root(), ts$annotation$psi))
psi_annot[, feature_id := psi_uid]

universe <- list()
for (mod in enabled) {
    spec <- ts$modalities[[mod]]
    if (mod == "expression_nearest_gene") {
        links <- window_links(gene_annot[, .(chrom, start, end, feature_id)],
                              spec$window_bp)
    } else if (mod == "psi") {
        links <- window_links(psi_annot[, .(chrom, start, end, feature_id)],
                              spec$window_bp)
    } else if (mod == "expression_abc") {
        ## ABC links are enhancer-to-gene calls from an external resource. The
        ## VMR side is resolved by overlap with the enhancer interval, so an ABC
        ## link exists only where a VMR sits in a called enhancer -- a much
        ## smaller and more specific universe than the distance windows.
        abc_f <- file.path(repo_root(), spec$link_source)
        if (!file.exists(abc_f)) {
            warning("ABC link source missing, skipping modality: ", abc_f)
            next
        }
        abc <- fread(abc_f)
        need <- c("PeakID_seqnames", "PeakID_start", "PeakID_end", "TargetGene")
        if (!all(need %in% names(abc))) {
            stop("ABC file lacks expected columns ",
                 paste(setdiff(need, names(abc)), collapse = ", "),
                 "; found: ", paste(head(names(abc), 12), collapse = ", "))
        }
        ## TargetGene is an Ensembl gene ID and matches the assay's gene_id once
        ## both are stripped of their version suffix. Verified 2026-09-02 that
        ## these coordinates are hg38: the file's TargetGeneTSS agrees with the
        ## gencode-v47 annotated TSS to a median of 1 bp (93.4% within 1 kb), so
        ## no liftover is applied.
        strip_ver <- function(x) sub("\\..*$", "", as.character(x))
        gene_annot[, gid := strip_ver(feature_id)]
        abc2 <- data.table(
            chrom = norm_chr(abc$PeakID_seqnames),
            start = as.integer(abc$PeakID_start),
            end   = as.integer(abc$PeakID_end),
            gid   = strip_ver(abc$TargetGene))
        abc2 <- merge(abc2, unique(gene_annot[, .(gid, feature_id)], by = "gid"),
                      by = "gid")
        if (nrow(abc2) == 0) stop("No ABC target gene matched the assay annotation")
        egr <- GRanges(abc2$chrom, IRanges(abc2$start, abc2$end))
        hits <- findOverlaps(vmr_gr, egr)
        links <- if (length(hits) == 0) data.table() else unique(data.table(
            vmr_id = mcols(vmr_gr)$vmr_id[queryHits(hits)],
            feature_id = abc2$feature_id[subjectHits(hits)],
            distance = 0L), by = c("vmr_id", "feature_id"))
    } else {
        stop("Unknown modality: ", mod)
    }

    if (nrow(links) == 0) {
        warning("Modality ", mod, " produced no links")
        next
    }
    links[, modality := mod]
    fwrite(links, file.path(run_dir, "links", paste0(mod, "-links.tsv.gz")),
           sep = "\t")
    universe[[mod]] <- data.table(
        modality = mod, label = spec$label,
        window_bp = spec$window_bp %||% NA_integer_,
        n_vmrs_linked = uniqueN(links$vmr_id),
        n_features_linked = uniqueN(links$feature_id),
        n_pairs = nrow(links),
        median_distance = median(links$distance),
        n_vmrs_in_catalog = nrow(vmr)
    )
    message("[07] ", mod, ": ", nrow(links), " pairs, ",
            uniqueN(links$vmr_id), " VMRs, ",
            uniqueN(links$feature_id), " features")
}

if (length(universe) == 0) stop("No modality produced any link")
uni <- rbindlist(universe, fill = TRUE)
uni[, `:=`(cohort = cohort, region = region, run_id = opts$run_id,
           vmr_set_id = mf("vmr_set_id", required = FALSE))]
write_atomic(uni, file.path(run_dir, "results", "tested-universe.tsv"))
print(uni)
