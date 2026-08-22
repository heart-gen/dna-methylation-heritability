#!/usr/bin/env Rscript

## Stage 05: fail-closed QC for the observed relative score. This gate assesses
## reconciliation, computational completeness, simulation-domain eligibility,
## and score integrity; it never re-evaluates or weakens the absolute-PVE gate.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(run_dir = ""))
if (!nzchar(cli$run_dir)) stop("--run-dir is required")
run_dir <- normalizePath(cli$run_dir)
manifest <- read_tsv(file.path(run_dir, "manifest.tsv"))
mval <- function(field) {
    value <- manifest$value[manifest$field == field]
    if (length(value) != 1L) stop("Run manifest lacks unique field: ", field)
    as.character(value[[1L]])
}
score_path <- file.path(
    run_dir, "results", "combined",
    sprintf("local-genetic-control-%s-%s-vmrs.tsv",
            mval("cohort"), mval("region"))
)
score <- read_tsv(score_path)
recon <- read_tsv(file.path(
    run_dir, "results", "combined", "task-reconciliation.tsv"
))
if (nrow(recon) != 1L) stop("Task reconciliation must contain one row")

parse_flag <- function(x, field) {
    value <- tolower(trimws(as.character(x)))
    if (any(is.na(x) | !value %in% c("true", "false", "t", "f", "1", "0"))) {
        stop("Invalid or missing ", field)
    }
    value %in% c("true", "t", "1")
}
feature_complete <- parse_flag(score$feature_complete, "feature_complete")
computational_failure <- parse_flag(
    score$computational_failure, "computational_failure"
)
eligible <- parse_flag(
    score$local_genetic_control_eligible, "local_genetic_control_eligible"
)
absolute_allowed <- parse_flag(
    score$absolute_pve_interpretation_allowed,
    "absolute_pve_interpretation_allowed"
)

threshold_lines <- readLines(file.path(run_dir, "config", "thresholds.yml"),
                             warn = FALSE)
outside_line <- grep("^[[:space:]]+max_outside_calibration_domain:",
                     threshold_lines, value = TRUE)
if (length(outside_line) != 1L) {
    stop("Cannot resolve max_outside_calibration_domain")
}
max_outside <- as_num(sub(".*:[[:space:]]*", "", outside_line),
                      "max_outside_calibration_domain")
complete_n <- sum(feature_complete & !computational_failure)
eligible_n <- sum(eligible)
eligible_rate <- if (complete_n) eligible_n / complete_n else 0
eligible_score <- score$local_snp_contribution_score[eligible]
raw_eligible <- score$pve_cis_joint_calibrated[eligible]
tie_counts <- table(raw_eligible, useNA = "no")
max_tie_fraction <- if (eligible_n && length(tie_counts)) {
    max(tie_counts) / eligible_n
} else NA_real_
lower_boundary <- upper_boundary <- rep(FALSE, nrow(score))
if (eligible_n) {
    lower_boundary[eligible] <- parse_flag(
        score$pve_lower_boundary_hit[eligible], "pve_lower_boundary_hit"
    )
    upper_boundary[eligible] <- parse_flag(
        score$pve_upper_boundary_hit[eligible], "pve_upper_boundary_hit"
    )
}
boundary <- eligible & (lower_boundary | upper_boundary)
boundary_rate <- if (eligible_n) sum(boundary) / eligible_n else NA_real_

checks <- data.frame(
    criterion = c(
        "task_reconciliation_complete", "zero_computational_failures",
        "within_domain_rate", "score_is_nondegenerate",
        "absolute_pve_is_prohibited", "interpretation_decision_is_locked"
    ),
    observed = c(
        as.numeric(recon$unique_task_rows == recon$expected &&
                   recon$unaccounted == 0 && recon$duplicate_task_ids == 0 &&
                   recon$unexpected_task_ids == 0),
        sum(computational_failure), eligible_rate,
        length(unique(eligible_score)), sum(absolute_allowed),
        as.numeric(all(score$local_genetic_control_decision ==
            "PASS_RELATIVE_GENETIC_CONTROL_FAIL_ABSOLUTE_LOCUS_PVE"))
    ),
    rule = c("equal", "equal", "greater_or_equal", "greater_or_equal",
             "equal", "equal"),
    threshold = c(1, 0, 1 - max_outside, 2, 0, 1),
    stringsAsFactors = FALSE
)
checks$pass <- with(checks,
    ifelse(rule == "equal", observed == threshold,
           observed >= threshold))
checks$detail <- c(
    sprintf("%d expected; %d unique; %d unaccounted",
            recon$expected, recon$unique_task_rows, recon$unaccounted),
    sprintf("%d rows flagged", sum(computational_failure)),
    sprintf("%d/%d complete-feature VMRs eligible", eligible_n, complete_n),
    sprintf("%d unique values; max tie fraction %.4f; boundary rate %.4f",
            length(unique(eligible_score)), max_tie_fraction, boundary_rate),
    "absolute_pve_interpretation_allowed must be FALSE for every row",
    paste(unique(score$local_genetic_control_decision), collapse = ",")
)
write_tsv(checks, file.path(
    run_dir, "results", "combined", "observed-score-qc.tsv"
))

all_pass <- all(checks$pass)
smoke_run <- identical(toupper(mval("smoke_run")), "TRUE")
decision <- if (!all_pass) {
    "FAIL_OBSERVED_RELATIVE_SCORE_QC"
} else if (smoke_run) {
    "PASS_SMOKE_ONLY_NOT_ACCEPTABLE"
} else {
    "PASS_RELATIVE_SCORE_OBSERVED_QC"
}
decision_table <- data.frame(
    run_id = mval("run_id"), cohort = mval("cohort"), region = mval("region"),
    vmr_set_id = mval("vmr_set_id"), decision = decision,
    relative_score_authorized = all_pass && !smoke_run,
    absolute_pve_authorized = FALSE,
    expected_vmrs = recon$expected, feature_complete_vmrs = complete_n,
    eligible_vmrs = eligible_n, eligible_rate = eligible_rate,
    max_tie_fraction = max_tie_fraction, boundary_rate = boundary_rate,
    stringsAsFactors = FALSE
)
write_tsv(decision_table, file.path(
    run_dir, "results", "combined", "observed-score-decision.tsv"
))
cat(decision, "\n")
if (!all_pass) quit(save = "no", status = 1L)
