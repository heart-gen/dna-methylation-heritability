#!/usr/bin/env Rscript
#### 03_local_snp_prediction -- pool out-of-fold predictions into metrics ####
##
## Usage:
##   Rscript _h/03_combine_oof.R --run-id lsp-AA-caudate-20260817
##
## Every metric in config/prediction.yml is computed from POOLED out-of-fold
## predictions. Two rules are enforced here rather than left to the analyst:
##
##  1. r2_pred_oof = 1 - SSE/SST, computed against the observed variance. It can
##     be negative, and negatives are RETAINED -- a locus whose predictor is
##     worse than the mean is a real and reportable result.
##  2. cor2_oof is computed and reported alongside, never substituted for
##     r2_pred_oof when r2_pred_oof is unfavorable (AGENTS.md 3). It is squared
##     correlation and is insensitive to scale and intercept error, which is
##     exactly why it flatters a miscalibrated predictor.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "03_local_snp_prediction"

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest_early <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval_early <- function(f) {
    v <- manifest_early$value[manifest_early$field == f]
    if (length(v) == 0) NA_character_ else v[1]
}

tasks <- fread(file.path(run_dir, "task-manifest.tsv"))
oof_dir <- file.path(run_dir, "results", "oof")
files <- list.files(oof_dir, pattern = "\\.tsv$", full.names = TRUE)
if (length(files) == 0) stop("No out-of-fold predictions under ", oof_dir)

preds <- rbindlist(lapply(files, fread), fill = TRUE)

## ------------------------------------------------------------ reconciliation
## AGENTS.md 9: zero tolerance for unexplained failures. Do this BEFORE any
## metric is computed, so a partial run cannot produce a plausible-looking table.
##
## Four outcomes, and every expected task must land in exactly one. A locus the
## upstream reader refused (no SNP in the cis window, fewer than min_cis_variants)
## is `qc_failed` with a recorded reason -- a documented exclusion, not a crash --
## and 02_fit_oof.R writes a `.skip` file saying which. Folding those into
## `failed` would make an ordinary genomic fact look like a computational defect;
## dropping them silently would leave the run unaccounted.
completed <- unique(as.character(preds$vmr_id))
skips <- list.files(oof_dir, pattern = "\\.skip$", full.names = TRUE)
skip_reasons <- if (length(skips)) {
    data.table(vmr_id = sub("\\.skip$", "", basename(skips)),
               detail = vapply(skips, function(f) paste(readLines(f), collapse = " "),
                               character(1), USE.NAMES = FALSE))
} else data.table(vmr_id = character(), detail = character())
qc_failed <- skip_reasons$vmr_id
failed <- sub("\\.txt$", "", list.files(file.path(run_dir, "results", "failures"),
                                        pattern = "\\.txt$"))

recon <- reconcile(expected = tasks$vmr_id, completed = completed,
                   qc_failed = qc_failed, failed = failed,
                   run = list(dir = run_dir),
                   allow_failures = identical(mval_early("smoke_run"), "TRUE"))
if (nrow(skip_reasons)) {
    write_atomic(skip_reasons,
                 file.path(run_dir, "results", "qc-failed-loci.tsv"))
}

comb_dir <- file.path(run_dir, "results", "combined")
dir.create(comb_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------------ metrics
metrics <- preds[, {
    sse <- sum((y_obs - y_pred)^2)
    sst <- sum((y_obs - mean(y_obs))^2)
    ## Calibration of a predictor: regress observed on predicted. Intercept 0
    ## and slope 1 is perfect calibration; slope < 1 is the usual overfitting
    ## signature (predictions too spread out).
    cal <- if (stats::sd(y_pred) > 0) {
        stats::coef(stats::lm(y_obs ~ y_pred))
    } else c(NA_real_, NA_real_)
    .(
        n_donors_predicted = uniqueN(donor),
        n_predictions      = .N,
        r2_pred_oof        = if (sst > 0) 1 - sse / sst else NA_real_,
        cor2_oof           = if (stats::sd(y_pred) > 0 && stats::sd(y_obs) > 0) {
                                 stats::cor(y_obs, y_pred)^2
                             } else NA_real_,
        rmse               = sqrt(mean((y_obs - y_pred)^2)),
        mae                = mean(abs(y_obs - y_pred)),
        calibration_intercept = cal[1],
        calibration_slope     = cal[2],
        screening_pass_frequency = mean(screened_in),
        median_n_variants  = stats::median(n_variants),
        median_alpha       = stats::median(alpha, na.rm = TRUE)
    )
}, by = vmr_id]

## Carry the upstream handles so a downstream module never has to guess which
## catalog these predictions belong to.
metrics <- merge(
    metrics,
    tasks[, .(vmr_id, chrom, start, end,
              local_genetic_control_eligible,
              local_genetic_control_exclusion_reason)],
    by = "vmr_id", all.x = TRUE
)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]; if (length(v) == 0) NA_character_ else v[1]
}
metrics[, `:=`(region = mval("region"), population = mval("cohort"),
               vmr_set_id = mval("vmr_set_id"),
               upstream_lgv_run_id = mval("upstream_local_genetic_variance_run_id"))]

## ------------------------------------------------- fold/repeat diagnostics
## AGENTS.md 7.3 requires fold and repeat diagnostics alongside the pooled
## metric. The headline r2_pred_oof pools every repeat into one SSE/SST, which
## is the right summary but hides the case that matters: a locus whose accuracy
## swings between repeats is being driven by a particular donor partition, not
## by a stable cis signal. Reported per repeat, and per fold within repeat.
r2_of <- function(obs, pred) {
    sst <- sum((obs - mean(obs))^2)
    if (sst > 0) 1 - sum((obs - pred)^2) / sst else NA_real_
}
by_repeat <- preds[, .(n_predictions = .N,
                       r2_pred_oof = r2_of(y_obs, y_pred),
                       screening_pass_frequency = mean(screened_in)),
                   by = .(vmr_id, repeat_i)]
by_fold <- preds[, .(n_predictions = .N,
                     r2_pred_oof = r2_of(y_obs, y_pred),
                     screened_in = as.logical(screened_in[1])),
                 by = .(vmr_id, repeat_i, outer_fold)]
write_atomic(by_repeat, file.path(comb_dir, "fold-repeat-diagnostics.tsv"))
write_atomic(by_fold, file.path(comb_dir, "fold-diagnostics.tsv"))

## Spread across repeats, carried onto the per-VMR table so the QC stage and any
## downstream consumer can see partition sensitivity without reopening the
## diagnostics.
spread <- by_repeat[, .(r2_repeat_min = min(r2_pred_oof, na.rm = TRUE),
                        r2_repeat_max = max(r2_pred_oof, na.rm = TRUE),
                        r2_repeat_sd = stats::sd(r2_pred_oof, na.rm = TRUE)),
                    by = vmr_id]
metrics <- merge(metrics, spread, by = "vmr_id", all.x = TRUE)

write_atomic(metrics, file.path(comb_dir,
    sprintf("oof-prediction-%s-%s-vmrs.tsv", mval("cohort"), mval("region"))))


## The BH diagnostic across VMRs. config/prediction.yml locks the screen at a
## fixed within-fold permutation alpha BECAUSE a fold cannot correct across a
## VMR set it cannot see; this is the reported view of how the screen behaved
## run-wide, never a gate that changed any prediction.
screen_tab <- preds[!is.na(screen_p),
                    .(min_screen_p = min(screen_p)), by = vmr_id]
if (nrow(screen_tab)) {
    screen_tab[, min_screen_p_bh := stats::p.adjust(min_screen_p, method = "BH")]
    write_atomic(screen_tab, file.path(comb_dir, "screen-bh-diagnostic.tsv"))
}

## Per-donor prediction counts: config/prediction.yml requires them, and an
## uneven count is the symptom of folds that silently dropped donors.
per_donor <- preds[, .(n_vmrs_predicted = uniqueN(vmr_id)), by = donor]
write_atomic(per_donor, file.path(comb_dir, "predictions-per-donor.tsv"))

## --------------------------------------------------------------------- QC
qc <- data.table(
    region = mval("region"), population = mval("cohort"),
    expected_vmrs = nrow(tasks),
    scored_vmrs = nrow(metrics),
    failed_vmrs = length(failed),
    median_r2_pred_oof = stats::median(metrics$r2_pred_oof, na.rm = TRUE),
    mean_r2_pred_oof = mean(metrics$r2_pred_oof, na.rm = TRUE),
    frac_r2_positive = mean(metrics$r2_pred_oof > 0, na.rm = TRUE),
    median_cor2_oof = stats::median(metrics$cor2_oof, na.rm = TRUE),
    median_calibration_slope = stats::median(metrics$calibration_slope, na.rm = TRUE),
    mean_screening_pass = mean(metrics$screening_pass_frequency, na.rm = TRUE),
    donor_prediction_counts_equal = uniqueN(per_donor$n_vmrs_predicted) == 1L,
    complete = length(failed) == 0L
)
write_atomic(qc, file.path(comb_dir, "prediction-run-qc.tsv"))
print(qc)

writeLines(capture.output(sessionInfo()), file.path(comb_dir, "session-info.txt"))

## A sanity flag, not a gate: if the median honest R2 comes back near the ~0.85
## the legacy code reported, something has leaked and the run must be audited
## before it is used.
if (isTRUE(qc$median_r2_pred_oof > 0.5)) {
    warning("median r2_pred_oof is ", round(qc$median_r2_pred_oof, 3),
            ", which is in the range defect E1 produced by leakage. ",
            "Audit the fold logic in 02_fit_oof.R before using this run.",
            call. = FALSE)
}
