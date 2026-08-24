#!/usr/bin/env Rscript
#### 04_repeat_repressive_architecture -- per-VMR annotation and covariate table ####
##
## Usage:
##   Rscript _h/01_build_features.R --run-id rra-AA-caudate-20260817
##
## Builds one row per interpretable VMR: the outcomes (H3K9me3 overlap,
## quiescent-chromatin overlap, LINE/L1 overlap), the primary predictor
## (standardized local SNP contribution score), the secondary predictor
## (r2_pred_oof), and
## every adjustment covariate named in config/repeat_annotations.yml.
##
## The covariate list is not decorative. VMR length, CpG count and density,
## mappability, segmental-duplication overlap and tested-SNP count all correlate
## with BOTH repeat content and the ability to estimate local genetic variance
## at all. An unadjusted overlap test recovers those confounds, not biology.
## This script therefore FAILS if a declared covariate cannot be constructed,
## rather than dropping it and proceeding with a thinner model.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(rtracklayer)
})

MODULE <- "04_repeat_repressive_architecture"

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]; if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mval("cohort"); region <- mval("region")
annot <- load_config("repeat_annotations")

## -------------------------------------------------------- primary predictor
lcg <- load_local_genetic_control(
    mval("upstream_local_genetic_variance_run_id"),
    region = region, cohort = cohort, eligible_only = TRUE
)

vmr <- GRanges(seqnames = lcg$chrom,
               ranges = IRanges(start = lcg$start, end = lcg$end),
               vmr_id = lcg$vmr_id)

feat <- data.table(
    vmr_id = lcg$vmr_id, chrom = lcg$chrom,
    start = lcg$start, end = lcg$end,
    local_snp_contribution_score = lcg$local_snp_contribution_score,
    local_snp_contribution_score_z = lcg$local_snp_contribution_score_z,
    local_snp_contribution_quartile = lcg$local_snp_contribution_quartile
)

## ------------------------------------------------------ secondary predictor
pred_run <- mval("upstream_local_snp_prediction_run_id")
if (!is.na(pred_run) && nzchar(pred_run)) {
    pf <- list.files(file.path(repo_root(), "03_local_snp_prediction", "_m",
                               "runs", pred_run, "results", "combined"),
                     pattern = "^oof-prediction-.*\\.tsv$", full.names = TRUE)
    if (length(pf) != 1) stop("Expected one OOF prediction table for ", pred_run)
    pr <- fread(pf[1])
    if ("r_squared_cv" %in% names(pr)) {
        stop("Prediction table carries the banned legacy metric r_squared_cv ",
             "(AGENTS.md 3)")
    }
    feat <- merge(feat, pr[, .(vmr_id, r2_pred_oof)], by = "vmr_id", all.x = TRUE)
    feat[, r2_pred_oof_z := as.numeric(scale(r2_pred_oof))]
}

## ------------------------------------------------------------------ outcomes
#' Fraction of a VMR covered by an annotation, and a binary any-overlap call.
#'
#' Both are kept: the fraction is the more informative outcome but is sensitive
#' to VMR length (hence the length covariate), while the binary call is what the
#' legacy analysis reported and is needed for the old-vs-new comparison.
overlap_features <- function(gr, bed_path, prefix) {
    if (is.null(bed_path)) {
        stop("config/repeat_annotations.yml has no track for '", prefix,
             "'. It is a PI-lock key and must be set before 04 runs.")
    }
    p <- if (startsWith(bed_path, "/")) bed_path else file.path(repo_root(), bed_path)
    if (!file.exists(p)) stop("Annotation track not found: ", p)
    ann <- reduce(import(p))
    hits <- findOverlaps(gr, ann)
    ## Intersect widths, summed per VMR, then divided by VMR width.
    inter <- pintersect(gr[queryHits(hits)], ann[subjectHits(hits)])
    cov <- tapply(width(inter), queryHits(hits), sum)
    frac <- numeric(length(gr))
    frac[as.integer(names(cov))] <- as.numeric(cov)
    frac <- pmin(frac / width(gr), 1)
    out <- data.table(frac, frac > 0)
    setnames(out, c(paste0(prefix, "_frac"), paste0(prefix, "_any")))
    out
}

feat <- cbind(feat,
    overlap_features(vmr, annot$chromatin$h3k9me3, "h3k9me3"),
    overlap_features(vmr, annot$chromatin$quiescent, "quiescent"),
    overlap_features(vmr, annot$repeatmasker$line_l1_bed, "line_l1"))

## ---------------------------------------------------------------- covariates
feat[, `:=`(
    vmr_length = end - start + 1L,
    cpg_count = lcg$n_cpgs,
    mean_methylation = lcg$mean_methylation,
    methylation_variance = lcg$methylation_variance,
    tested_snp_count = lcg$n_variants
)]
feat[, cpg_density := cpg_count / vmr_length]

## Mappability and segdup overlap come from the same sensitivity assets, because
## they must be identical between the primary model and the restriction.
feat <- cbind(feat, overlap_features(vmr, annot$sensitivities$exclude_segdups$bed,
                                     "segdup"))

umap <- annot$sensitivities$high_mappability$umap_k24_bw
if (is.null(umap)) stop("No mappability track configured; it is a required covariate")
umap_p <- if (startsWith(umap, "/")) umap else file.path(repo_root(), umap)
if (!file.exists(umap_p)) stop("Mappability bigWig not found: ", umap_p)
feat[, mappability := {
    sc <- summary(BigWigFile(umap_p), which = vmr, type = "mean")
    as.numeric(unlist(lapply(sc, function(x) if (length(x$score)) x$score else NA_real_)))
}]

## Fail loudly on anything the declared model needs and we could not build.
required <- c("vmr_length", "cpg_count", "cpg_density", "mean_methylation",
              "methylation_variance", "tested_snp_count", "mappability",
              "segdup_frac")
missing_cols <- setdiff(required, names(feat))
if (length(missing_cols) > 0) {
    stop("Declared covariates could not be constructed: ",
         paste(missing_cols, collapse = ", "),
         "\n  Fix the upstream columns rather than fitting a thinner model.")
}
incomplete <- required[vapply(required, function(k) all(is.na(feat[[k]])), logical(1))]
if (length(incomplete) > 0) {
    stop("Covariate(s) constructed but entirely NA: ",
         paste(incomplete, collapse = ", "))
}

feat[, `:=`(region = region, population = cohort, vmr_set_id = mval("vmr_set_id"))]
write_atomic(feat, file.path(run_dir, "results", "vmr-features.tsv"))
message("[04] built features for ", nrow(feat), " interpretable VMRs")
