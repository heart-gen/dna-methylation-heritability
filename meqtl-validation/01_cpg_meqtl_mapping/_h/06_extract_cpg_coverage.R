#!/usr/bin/env Rscript
# Extract real per-CpG coverage from staff-repo stats.rda BSobj for prepared VMR CpGs.
#
# Usage:
#   Rscript 06_extract_cpg_coverage.R --region caudate --chrom 22
#
# Writes cpg_coverage.chrN.tsv and updates cpg_vmr_map.chrN.tsv with coverage columns.

suppressPackageStartupMessages({
  library(bsseq)
  library(DelayedMatrixStats)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit[1] == length(args)) stop("Missing value for ", flag)
  args[hit[1] + 1]
}

region <- get_arg("--region")
chrom_raw <- get_arg("--chrom")
min_cov <- as.integer(get_arg("--min-coverage", "5"))
project_root <- get_arg("--project-root",
  "/projects/b1213/users/kynon/projects/dna-methylation-heritability")
staff_root <- get_arg("--staff-root",
  "/projects/b1213/users/alexis/projects/dna-methylation-heritability")
if (is.null(region) || is.null(chrom_raw)) {
  stop("Usage: Rscript 06_extract_cpg_coverage.R --region REGION --chrom N [--min-coverage 5]")
}

chrom <- gsub("^chr", "", chrom_raw)
chrom_label <- paste0("chr", chrom)

prep_dir <- file.path(project_root, "meqtl-validation", "01_cpg_meqtl_mapping",
                      region, "_m", "prepared")
preflight <- file.path(project_root, "meqtl-validation", "01_cpg_meqtl_mapping",
                       region, "_m", "preflight", "sample_inclusion_primary.tsv")
stats_rda <- file.path(staff_root, "vmr-analysis", region, "_m", "cpg",
                       paste0("chr_", chrom), "stats.rda")
map_path <- file.path(prep_dir, paste0("cpg_vmr_map.", chrom_label, ".tsv"))
out_path <- file.path(prep_dir, paste0("cpg_coverage.", chrom_label, ".tsv"))

if (!file.exists(stats_rda)) stop("Missing stats.rda: ", stats_rda)
if (!file.exists(preflight)) stop("Missing preflight inclusion list: ", preflight)
dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)

inc <- fread(preflight)
keep_br <- unique(inc$brnum)

message("Loading ", stats_rda)
load(stats_rda)  # BSobj, means, sds
if (!exists("BSobj")) stop("BSobj not found in stats.rda")

br <- as.character(colData(BSobj)$brnum)
keep_idx <- which(br %in% keep_br)
if (length(keep_idx) == 0) {
  message("WARNING: no overlap with inclusion list; using all BSobj samples")
  keep_idx <- seq_len(ncol(BSobj))
}
BSobj <- BSobj[, keep_idx]
n_samp <- ncol(BSobj)

cov <- getCoverage(BSobj, type = "Cov")
pos <- start(BSobj)

if (file.exists(map_path)) {
  cmap <- fread(map_path)
  target_pos <- as.integer(cmap$pos_1based)
  message("Annotating ", length(target_pos), " prepared VMR CpGs")
} else {
  target_pos <- as.integer(pos)
  cmap <- data.table(
    phenotype_id = paste0(chrom_label, "_", target_pos),
    chrom = chrom_label,
    pos_1based = target_pos
  )
  message("No cpg_vmr_map; exporting coverage for all ", length(target_pos), " loci")
}

idx <- match(target_pos, pos)
found <- which(!is.na(idx))
miss <- which(is.na(idx))
message("matched ", length(found), " / ", length(target_pos),
        " (unmatched ", length(miss), ")")

mean_cov <- rep(NA_real_, length(target_pos))
median_cov <- rep(NA_real_, length(target_pos))
n_ge <- rep(NA_integer_, length(target_pos))
frac_ge <- rep(NA_real_, length(target_pos))

if (length(found) > 0) {
  rows <- idx[found]
  chunk <- 50000L
  for (start in seq(1L, length(rows), by = chunk)) {
    end <- min(length(rows), start + chunk - 1L)
    rr <- rows[start:end]
    sub <- as.matrix(cov[rr, , drop = FALSE])
    mean_cov[found[start:end]] <- rowMeans(sub)
    median_cov[found[start:end]] <- matrixStats::rowMedians(sub)
    n_ge[found[start:end]] <- as.integer(rowSums(sub >= min_cov))
    frac_ge[found[start:end]] <- n_ge[found[start:end]] / n_samp
    message("  coverage chunk ", start, "-", end)
  }
}

out <- data.table(
  phenotype_id = paste0(chrom_label, "_", target_pos),
  chrom = chrom_label,
  pos_1based = target_pos,
  n_samples_coverage = n_samp,
  min_coverage_threshold = min_cov,
  mean_coverage = mean_cov,
  median_coverage = median_cov,
  n_samples_cov_ge_min = n_ge,
  fraction_samples_cov_ge_min = frac_ge,
  coverage_matched_in_bsobj = !is.na(idx)
)
fwrite(out, out_path, sep = "\t")
message("Wrote ", out_path)

if (file.exists(map_path)) {
  cmap2 <- copy(cmap)
  drop_cols <- intersect(
    names(cmap2),
    c("n_samples_coverage", "min_coverage_threshold", "mean_coverage",
      "median_coverage", "n_samples_cov_ge_min", "fraction_samples_cov_ge_min",
      "coverage_matched_in_bsobj", "coverage_source")
  )
  if (length(drop_cols)) cmap2[, (drop_cols) := NULL]
  merged <- merge(
    cmap2, out,
    by = c("phenotype_id", "chrom", "pos_1based"),
    all.x = TRUE, sort = FALSE
  )
  merged[, coverage_source := "stats.rda_BSobj_getCoverage"]
  fwrite(merged, map_path, sep = "\t")
  message("Updated ", map_path)
}

sum_path <- file.path(prep_dir, paste0("cpg_coverage_summary.", chrom_label, ".tsv"))
fwrite(data.table(
  region = region,
  chrom = chrom_label,
  n_target_cpgs = length(target_pos),
  n_matched = length(found),
  n_unmatched = length(miss),
  n_samples_coverage = n_samp,
  min_coverage_threshold = min_cov,
  mean_of_mean_coverage = mean(mean_cov, na.rm = TRUE),
  median_of_mean_coverage = median(mean_cov, na.rm = TRUE),
  mean_frac_samples_cov_ge_min = mean(frac_ge, na.rm = TRUE),
  out_path = out_path
), sum_path, sep = "\t")
message("Wrote ", sum_path)
