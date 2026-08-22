#!/usr/bin/env Rscript

## Observed-regime grid, stage C: reconcile, apply the frozen joint model, and
## answer the two questions the AR(1) validation grid could not.
##
##   Q1  Is the observed 62.7% lower-boundary rate reproduced once the
##       estimator is run on real cis-window genotypes at the observed n?
##   Q2  Inside the lower-boundary set, does the unbounded estimate still carry
##       rank information about true local PVE?
##
## The frozen model is applied unmodified; nothing here refits anything.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
h_dir <- dirname(script_path)
source(file.path(h_dir, "00_functions.R"))
source(file.path(h_dir, "joint_pve_functions.R"))

cli <- parse_cli(list(run_dir = "", bootstrap = "2000"))
if (!nzchar(cli$run_dir)) stop("--run-dir is required")
run_dir <- normalizePath(cli$run_dir)
n_boot <- as_int(cli$bootstrap, "bootstrap")
combined_dir <- file.path(run_dir, "results", "combined")
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read_tsv(file.path(run_dir, "manifest.tsv"))
mval <- function(field) {
    value <- manifest$value[manifest$field == field]
    if (length(value) != 1L) stop("Run manifest lacks unique field: ", field)
    as.character(value[[1L]])
}
scenarios <- read_tsv(file.path(run_dir, "config", "scenario-manifest.tsv"))
row_dir <- file.path(run_dir, "results", "scenario_rows")
row_files <- list.files(row_dir, pattern = "^scenario-[0-9]{7}\\.tsv$",
                        full.names = TRUE)
if (!length(row_files)) stop("No scenario rows to combine")
rows <- do.call(rbind, lapply(row_files, read_tsv))
rows <- rows[order(rows$scenario_id), , drop = FALSE]
write_tsv(rows, file.path(combined_dir, "observed-regime-features.tsv"))

## ---- reconciliation -------------------------------------------------------
expected <- as.integer(scenarios$scenario_id)
observed_ids <- as.integer(rows$scenario_id)
reconciliation <- data.frame(
    expected_scenarios = length(expected),
    scenario_files = length(row_files),
    unique_scenario_rows = length(unique(observed_ids)),
    completed = sum(rows$terminal_status %in% "completed"),
    qc_failed = sum(rows$terminal_status %in% "qc_failed"),
    excluded = sum(rows$terminal_status %in% "excluded"),
    computational_failure = sum(rows$computational_failure %in% TRUE),
    duplicate = sum(duplicated(observed_ids)),
    unexpected = length(setdiff(observed_ids, expected)),
    unaccounted = length(setdiff(expected, observed_ids)),
    stringsAsFactors = FALSE
)
write_tsv(reconciliation, file.path(combined_dir, "regime-reconciliation.tsv"))

## ---- apply the frozen joint model ----------------------------------------
model_path <- mval("joint_model_path")
observed_sha <- tolower(sub(" .*$", "", system2(
    "sha256sum", normalizePath(model_path), stdout = TRUE)[[1L]]))
if (!identical(observed_sha, tolower(mval("joint_model_sha256")))) {
    stop("Frozen joint-model checksum changed after stage A")
}
model <- readRDS(model_path)
if (!identical(model$gate_version, "final_joint_pve_v1")) {
    stop("Unexpected joint-model gate version: ", model$gate_version)
}
complete <- which(rows$feature_complete %in% TRUE &
                  !(rows$computational_failure %in% TRUE))
if (!length(complete)) stop("No complete observed-regime feature rows")
predicted <- predict_joint_pve(model, rows[complete, , drop = FALSE])
if (nrow(predicted) != length(complete)) stop("Prediction row-count mismatch")
for (field in names(predicted)) rows[[field]] <- NA
for (field in names(predicted)) rows[[field]][complete] <- predicted[[field]]
rows$joint_model_run_id <- mval("joint_model_run_id")
rows$joint_model_sha256 <- observed_sha
rows$absolute_pve_interpretation_allowed <- FALSE
write_tsv(rows, file.path(combined_dir, "observed-regime-estimates.tsv"))

fit <- rows[complete, , drop = FALSE]
if (!all(c("pve_cis_joint_unbounded", "pve_cis_joint_calibrated",
           "pve_lower_boundary_hit") %in% names(fit))) {
    stop("Frozen model did not emit the expected estimate fields")
}
unbounded <- as.numeric(fit$pve_cis_joint_unbounded)
calibrated <- as.numeric(fit$pve_cis_joint_calibrated)
truth <- as.numeric(fit$true_h2)
lower_hit <- fit$pve_lower_boundary_hit %in% TRUE

spearman <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L || length(unique(x[ok])) < 2L ||
        length(unique(y[ok])) < 2L) return(NA_real_)
    suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
}
kendall <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L || length(unique(x[ok])) < 2L ||
        length(unique(y[ok])) < 2L) return(NA_real_)
    suppressWarnings(stats::cor(x[ok], y[ok], method = "kendall"))
}

## ---- Q1: boundary rate, overall and stratified ---------------------------
boundary_row <- function(label, level, index) {
    data.frame(
        stratum = label, level = as.character(level),
        n_scenarios = sum(index),
        lower_boundary_rate = mean(lower_hit[index]),
        median_unbounded = stats::median(unbounded[index]),
        min_unbounded = min(unbounded[index]),
        stringsAsFactors = FALSE
    )
}
boundary_tables <- list(boundary_row("overall", "all", rep(TRUE, nrow(fit))))
for (field in c("true_h2", "architecture", "num_snps_stratum",
                "p_eff_stratum")) {
    for (level in sort(unique(fit[[field]]))) {
        boundary_tables[[length(boundary_tables) + 1L]] <-
            boundary_row(field, level, fit[[field]] == level)
    }
}
## Low-PVE cells are the ones comparable to the bulk of real VMRs.
low <- truth <= 0.1
boundary_tables[[length(boundary_tables) + 1L]] <-
    boundary_row("true_h2_at_most_0.1", "pooled", low)
boundary <- do.call(rbind, boundary_tables)
write_tsv(boundary, file.path(combined_dir, "regime-boundary-rate.tsv"))

## ---- Q2: rank information overall and below the boundary -----------------
sub_index <- which(lower_hit)
sub_truth <- truth[sub_index]
sub_unbounded <- unbounded[sub_index]
boot <- rep(NA_real_, n_boot)
if (length(sub_index) >= 10L && length(unique(sub_truth)) > 1L) {
    set.seed(as_int(mval("base_seed"), "base_seed") + 7L)
    for (b in seq_len(n_boot)) {
        pick <- sample(seq_along(sub_index), length(sub_index), replace = TRUE)
        boot[[b]] <- spearman(sub_unbounded[pick], sub_truth[pick])
    }
}
auc <- NA_real_
wilcox_p <- NA_real_
if (length(sub_index) && any(sub_truth == 0) && any(sub_truth > 0)) {
    a <- sub_unbounded[sub_truth > 0]
    b <- sub_unbounded[sub_truth == 0]
    auc <- mean(outer(a, b, ">")) + 0.5 * mean(outer(a, b, "=="))
    wilcox_p <- suppressWarnings(stats::wilcox.test(a, b)$p.value)
}
rank_info <- data.frame(
    metric = c(
        "spearman_overall_unbounded", "spearman_overall_calibrated",
        "spearman_low_h2_unbounded", "spearman_low_h2_calibrated",
        "n_lower_boundary", "unique_unbounded_below_boundary",
        "spearman_below_boundary", "kendall_below_boundary",
        "bootstrap_ci_lower", "bootstrap_ci_upper",
        "auc_nonnull_over_null_below_boundary", "wilcoxon_p_below_boundary"
    ),
    value = c(
        spearman(unbounded, truth), spearman(calibrated, truth),
        spearman(unbounded[low], truth[low]),
        spearman(calibrated[low], truth[low]),
        length(sub_index), length(unique(sub_unbounded)),
        spearman(sub_unbounded, sub_truth), kendall(sub_unbounded, sub_truth),
        stats::quantile(boot, 0.025, na.rm = TRUE),
        stats::quantile(boot, 0.975, na.rm = TRUE),
        auc, wilcox_p
    ),
    stringsAsFactors = FALSE
)
write_tsv(rank_info, file.path(combined_dir, "regime-rank-information.tsv"))

## ---- verdicts -------------------------------------------------------------
observed_boundary_rate <- NA_real_
observed_decision <- file.path(
    dirname(run_dir), mval("observed_run_id"),
    "results", "combined", "observed-score-decision.tsv"
)
if (file.exists(observed_decision)) {
    decision <- read_tsv(observed_decision)
    if ("boundary_rate" %in% names(decision)) {
        observed_boundary_rate <- as.numeric(decision$boundary_rate[[1L]])
    }
}
low_rate <- mean(lower_hit[low])
ci_lower <- rank_info$value[rank_info$metric == "bootstrap_ci_lower"]
ci_upper <- rank_info$value[rank_info$metric == "bootstrap_ci_upper"]
sub_rank_informative <- isTRUE(is.finite(ci_lower) && is.finite(ci_upper) &&
                                   ci_lower > 0)
rate_reproduced <- isTRUE(is.finite(observed_boundary_rate) &&
                              observed_boundary_rate <= max(boundary$lower_boundary_rate))
verdict <- data.frame(
    field = c(
        "run_id", "n_loci", "n_scenarios_complete",
        "observed_boundary_rate", "regime_boundary_rate_overall",
        "regime_boundary_rate_low_h2", "regime_boundary_rate_max_cell",
        "boundary_rate_reproduced", "sub_boundary_rank_informative",
        "unaccounted_scenarios", "computational_failures"
    ),
    value = as.character(c(
        mval("run_id"), mval("n_loci"), nrow(fit),
        observed_boundary_rate, mean(lower_hit), low_rate,
        max(boundary$lower_boundary_rate), rate_reproduced,
        sub_rank_informative, reconciliation$unaccounted,
        reconciliation$computational_failure
    )),
    stringsAsFactors = FALSE
)
write_tsv(verdict, file.path(combined_dir, "regime-verdict.tsv"))

cat("Observed-regime grid summary for", mval("run_id"), "\n")
cat("  complete scenarios          :", nrow(fit), "\n")
cat("  unaccounted scenarios       :", reconciliation$unaccounted, "\n")
cat("  observed boundary rate      :", observed_boundary_rate, "\n")
cat("  regime boundary rate        :", mean(lower_hit), "\n")
cat("  regime boundary rate h2<=0.1:", low_rate, "\n")
cat("  boundary rate reproduced    :", rate_reproduced, "\n")
cat("  sub-boundary rank info      :", sub_rank_informative, "\n")
