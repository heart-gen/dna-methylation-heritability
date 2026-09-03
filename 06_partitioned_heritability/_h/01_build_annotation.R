#!/usr/bin/env Rscript
#### 06_partitioned_heritability -- build the continuous annotation BED ####
##
## Usage:
##   Rscript _h/01_build_annotation.R --run-id sldsc-AA-caudate-YYYYMMDD
##
## Emits an hg38 BED whose 4th column is the continuous local SNP contribution
## z-score. 02_liftover_annotation.py lifts it to hg19 and 03_make_annot.py maps
## it onto the LD reference SNPs.
##
## This is the script that enforces AGENTS.md 3. The retired v1 analysis
## (region_heritability.py) did three banned things in eight lines: it filtered
## on `r_squared_cv >= 0.3`, it carried `h2_unscaled` as the annotation value,
## and it cut that value into quintile BEDs whose 4th column was then dropped
## entirely -- so the "continuous" annotation was neither continuous nor
## derived from an admissible quantity. None of that can be reached from here:
## the loader below refuses banned columns, and the checks refuse a
## thresholded, grouped, or degenerate annotation.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE", file.path(Sys.getenv("V2_REPO_ROOT", "."), "06_partitioned_heritability", "_h")), "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "06_partitioned_heritability"

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
## Select OUTSIDE the data.table `[`. Inside it the bare name `field` resolves
## to the COLUMN rather than to this argument, so `manifest$field == field` is
## column == column -- TRUE for every row -- and mf() silently returns the first
## row of the manifest whatever you asked for. Same trap as the one documented
## in 00_shared/gates.R::require_accepted_upstream().
mf <- function(field, required = TRUE) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0 || is.na(v[1])) {
        if (required) stop("Manifest has no value for '", field, "'")
        return(NA_character_)
    }
    v[1]
}

cohort <- mf("cohort")
region <- mf("region")
smoke  <- identical(mf("smoke_run"), "TRUE")
upstream_run <- mf("upstream_local_genetic_variance_run_id", required = !smoke)

ph <- load_run_config("partitioned_heritability", run_dir)
score_col <- ph$annotation$score_column

if (is.na(upstream_run)) {
    stop("This run has no upstream Module 02 run ID recorded. A smoke run may ",
         "skip the acceptance gate but still needs a score table to read.")
}

## load_local_genetic_control() is the shared loader used by Modules 03-05. It
## already refuses a table carrying h2_unscaled / r_squared_cv / h2_en_calibrated,
## refuses one that authorises absolute-PVE interpretation, and refuses one whose
## terminal decision is not the relative-score decision. Reusing it here means
## Module 06 cannot drift away from the interpretation rules the other modules
## follow.
dt <- load_local_genetic_control(
    upstream_run_id = upstream_run, region = region, cohort = cohort,
    eligible_only = isTRUE(ph$annotation$eligible_only)
)

if (!score_col %in% names(dt)) {
    stop("Score column '", score_col, "' absent from the Module 02 table")
}

## Belt and braces over the loader: the config names the banned columns
## explicitly so a future edit to either list is caught by the other.
banned <- intersect(ph$annotation$forbidden_columns, names(dt))
if (length(banned)) {
    stop("Module 02 table carries banned column(s): ",
         paste(banned, collapse = ", "), " (AGENTS.md 3)")
}

bed <- dt[, .(chrom = as.character(chrom),
              start = as.integer(start),
              end   = as.integer(end),
              score = as.numeric(get(score_col)),
              vmr_id = as.character(vmr_id))]

bed <- bed[is.finite(score)]
if (nrow(bed) == 0) stop("No VMRs with a finite ", score_col)

## LDSC's baseline model is autosomal, and a sex-chromosome VMR has no LD score
## to be partitioned against. Excluding them is a policy exclusion, recorded as
## such rather than silently dropped.
autosomes <- paste0("chr", 1:22)
bed[, chrom := ifelse(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))]
n_before <- nrow(bed)
excluded <- bed[!chrom %in% autosomes]
bed <- bed[chrom %in% autosomes]
if (nrow(excluded) > 0) {
    write_atomic(excluded, file.path(run_dir, "excluded",
                                     "non-autosomal-vmrs.tsv"))
}

if (any(bed$start >= bed$end)) stop("Degenerate VMR interval(s) in the score table")

## --- The annotation must be continuous -------------------------------------
## A thresholded or grouped annotation shows up as a tiny number of distinct
## values. The v1 quintile form would land at exactly 5. Refuse anything that
## looks like a partition rather than a gradient.
n_distinct <- uniqueN(bed$score)
if (n_distinct < 100) {
    stop("Annotation takes only ", n_distinct, " distinct values, which is a ",
         "partition and not a continuous score. AGENTS.md 3 bans the ",
         "thresholded v1 form of this annotation.")
}
sd_score <- sd(bed$score)
if (!is.finite(sd_score) || sd_score < ph$gates$min_annotation_sd) {
    stop("Annotation has collapsed to a constant (sd = ", sd_score, ")")
}
min_vmrs <- if (smoke) 1L else as.integer(ph$annotation$min_vmrs)
if (nrow(bed) < min_vmrs) {
    stop("Only ", nrow(bed), " autosomal VMRs carry a score; ",
         "config requires at least ", min_vmrs)
}

setorder(bed, chrom, start)
write_atomic(bed[, .(chrom, start, end, score)],
             file.path(run_dir, "annotation", "annotation-hg38.bed"),
             col.names = FALSE)
write_atomic(bed, file.path(run_dir, "annotation", "annotation-hg38-keyed.tsv"))

summary_dt <- data.table(
    cohort = cohort, region = region,
    vmr_set_id = mf("vmr_set_id", required = FALSE),
    upstream_local_genetic_variance_run_id = upstream_run,
    score_column = score_col,
    n_vmrs_scored = n_before,
    n_vmrs_autosomal = nrow(bed),
    n_vmrs_excluded_non_autosomal = nrow(excluded),
    n_distinct_values = n_distinct,
    score_mean = mean(bed$score), score_sd = sd_score,
    score_min = min(bed$score), score_max = max(bed$score),
    annotation_continuous = TRUE,
    annotation_thresholded = FALSE,
    absolute_pve_used = FALSE
)
write_atomic(summary_dt, file.path(run_dir, "annotation", "annotation-summary.tsv"))

message("[06] annotation: ", nrow(bed), " autosomal VMRs, ",
        n_distinct, " distinct values, sd ", signif(sd_score, 4))
