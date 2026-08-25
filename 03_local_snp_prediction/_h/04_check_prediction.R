#!/usr/bin/env Rscript
#### 03_local_snp_prediction -- acceptance gate ####
##
## Usage:
##   Rscript _h/04_check_prediction.R --run-id lsp-AA-caudate-20260823
##
## Mirrors Module 02 Stage 05. A completed SLURM job is not acceptance; this is.
## Every criterion is evaluated, all of them are written out, and the decision is
## the AND of them -- a gate that stops at the first failure hides the rest.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

MODULE <- "03_local_snp_prediction"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]; if (length(v) == 0) NA_character_ else v[1]
}
smoke <- identical(mval("smoke_run"), "TRUE")

comb_dir <- file.path(run_dir, "results", "combined")
mf <- list.files(comb_dir, pattern = "^oof-prediction-.*-vmrs\\.tsv$", full.names = TRUE)
if (length(mf) != 1) stop("Expected exactly one OOF prediction table in ", comb_dir)
metrics <- fread(mf)
qc <- fread(file.path(comb_dir, "prediction-run-qc.tsv"))
recon <- fread(file.path(run_dir, "task_reconciliation.tsv"))
rn <- function(k) as.integer(recon$n[recon$category == k])

prediction <- load_config("prediction")
tripwire <- tryCatch(config_get(prediction, "leakage_tripwire.max_median_r2_pred_oof"),
                     error = function(e) 0.5)

## ---------------------------------------------------------------- criteria
## The banned columns are the point of the module. r_squared_cv is defect E1's
## in-sample statistic and h2_en_calibrated is Module 02's retired estimator;
## neither may travel downstream under any name.
banned <- intersect(c("r_squared_cv", "h2_en_calibrated", "h2_unscaled"),
                    names(metrics))

per_donor <- fread(file.path(comb_dir, "predictions-per-donor.tsv"))

criteria <- data.table(
    criterion = c(
        "reconciliation_complete",
        "zero_computational_failures",
        "no_banned_legacy_columns",
        "leakage_tripwire_not_tripped",
        "finite_nondegenerate_r2",
        "negative_r2_retained",
        "donor_prediction_counts_equal",
        "screen_frequency_recorded"
    ),
    passed = c(
        rn("unaccounted") == 0L && rn("unexpected") == 0L,
        rn("failed") == 0L,
        length(banned) == 0L,
        !is.finite(qc$median_r2_pred_oof) ||
            qc$median_r2_pred_oof <= tripwire,
        sum(is.finite(metrics$r2_pred_oof)) > 0 &&
            stats::sd(metrics$r2_pred_oof, na.rm = TRUE) > 0,
        ## AGENTS.md 7.3: negative r2_pred_oof must be RETAINED, never clipped
        ## and never swapped for cor2_oof. If nothing is negative in a real run
        ## that is suspicious, but the testable property is that no clipping
        ## floor was applied -- no pile-up exactly at zero.
        sum(metrics$r2_pred_oof == 0, na.rm = TRUE) <= 1L,
        isTRUE(as.logical(qc$donor_prediction_counts_equal)),
        all(is.finite(metrics$screening_pass_frequency))
    ),
    detail = c(
        sprintf("unaccounted=%d unexpected=%d", rn("unaccounted"), rn("unexpected")),
        sprintf("failed=%d qc_failed=%d", rn("failed"), rn("qc_failed")),
        if (length(banned)) paste(banned, collapse = ",") else "none",
        sprintf("median_r2_pred_oof=%.4f threshold=%.2f",
                qc$median_r2_pred_oof, tripwire),
        sprintf("n_finite=%d sd=%.4f", sum(is.finite(metrics$r2_pred_oof)),
                stats::sd(metrics$r2_pred_oof, na.rm = TRUE)),
        sprintf("n_negative=%d n_exactly_zero=%d",
                sum(metrics$r2_pred_oof < 0, na.rm = TRUE),
                sum(metrics$r2_pred_oof == 0, na.rm = TRUE)),
        sprintf("distinct_counts=%d", uniqueN(per_donor$n_vmrs_predicted)),
        sprintf("mean_screening_pass=%.4f", qc$mean_screening_pass)
    )
)

decision <- if (all(criteria$passed)) {
    if (smoke) "PASS_SMOKE_ONLY_NOT_ACCEPTABLE" else "PASS_OOF_PREDICTION_QC"
} else {
    "FAIL_OOF_PREDICTION_QC"
}

criteria[, `:=`(run_id = opts$run_id, region = mval("region"),
                population = mval("cohort"), decision = decision)]
write_atomic(criteria, file.path(comb_dir, "prediction-qc-criteria.tsv"))
write_atomic(data.table(run_id = opts$run_id, decision = decision,
                        smoke_run = smoke,
                        n_criteria = nrow(criteria),
                        n_passed = sum(criteria$passed)),
             file.path(comb_dir, "prediction-decision.tsv"))
print(criteria[, .(criterion, passed, detail)])
message("[03] decision: ", decision)

if (!all(criteria$passed)) {
    quit(status = 1L)
}
