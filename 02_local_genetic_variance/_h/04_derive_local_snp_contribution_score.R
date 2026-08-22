#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(run_dir = "", input = "", output = ""))
if (nzchar(cli$run_dir)) {
    run_dir <- normalizePath(cli$run_dir)
    manifest <- read_tsv(file.path(run_dir, "manifest.tsv"))
    mval <- function(field) {
        value <- manifest$value[manifest$field == field]
        if (length(value) != 1L) stop("Run manifest lacks unique field: ", field)
        as.character(value[[1L]])
    }
    cli$input <- file.path(
        run_dir, "results", "combined", "observed-joint-estimates.tsv"
    )
    cli$output <- file.path(
        run_dir, "results", "combined",
        sprintf("local-genetic-control-%s-%s-vmrs.tsv",
                mval("cohort"), mval("region"))
    )
}
if (!nzchar(cli$input) || !nzchar(cli$output)) {
    stop("--run-dir or both --input and --output are required")
}

d <- read_tsv(cli$input)
banned <- intersect(c("h2_unscaled", "r_squared_cv"), names(d))
if (length(banned)) {
    stop("Banned legacy columns in score input: ", paste(banned, collapse = ", "))
}
required <- c(
    "vmr_id", "cohort", "region", "vmr_set_id",
    "chrom", "start", "end", "n_cpgs", "n_variants",
    "mean_methylation", "methylation_variance",
    "pve_cis_joint_unbounded", "pve_cis_joint_calibrated", "feature_complete",
    "computational_failure", "joint_pve_domain_status"
)
missing <- setdiff(required, names(d))
if (length(missing)) stop("Score input is missing: ", paste(missing, collapse = ", "))
if (anyDuplicated(d$vmr_id)) stop("Duplicate vmr_id in score input")
if (length(unique(d$cohort)) != 1L || length(unique(d$region)) != 1L ||
    length(unique(d$vmr_set_id)) != 1L) {
    stop("Score input must describe exactly one cohort, region, and vmr_set_id")
}

as_flag <- function(x, field) {
    value <- tolower(trimws(as.character(x)))
    if (any(is.na(x) | !value %in% c("true", "false", "t", "f", "1", "0"))) {
        stop("Invalid or missing ", field, " value")
    }
    value %in% c("true", "t", "1")
}
feature_complete <- as_flag(d$feature_complete, "feature_complete")
computational_failure <- as_flag(d$computational_failure, "computational_failure")
## The relative score ranks on the pre-clip estimate. Clipping is right for a
## descriptive PVE -- negative variance explained is not interpretable -- but
## it destroys ordering across the lower-boundary mass, which the 2026-08-22
## observed-regime grid showed is not noise: below the boundary the unbounded
## estimate tracks true PVE at Spearman 0.584 (bootstrap 95% CI 0.563-0.603,
## AUC 0.741). Clipping costs about 80% of the recoverable low-PVE ordering
## (Spearman 0.312 unbounded versus 0.059 clipped at true h2 <= 0.1).
## pve_cis_joint_calibrated is retained unchanged as the bounded descriptive
## output; neither column authorises absolute-PVE language.
score_basis <- "pve_cis_joint_unbounded"
finite_estimate <- is.finite(as.numeric(d[[score_basis]]))
domain_status <- as.character(d$joint_pve_domain_status)
within_domain <- !is.na(domain_status) & domain_status == "within_domain"
eligible <- feature_complete & !computational_failure & finite_estimate & within_domain

reason <- rep(NA_character_, nrow(d))
reason[!feature_complete] <- "incomplete_joint_features"
reason[feature_complete & computational_failure] <- "computational_failure"
reason[feature_complete & !computational_failure & !finite_estimate] <-
    "nonfinite_joint_estimate"
reason[feature_complete & !computational_failure & finite_estimate &
       is.na(domain_status)] <- "missing_joint_pve_domain_status"
reason[feature_complete & !computational_failure & finite_estimate &
       !is.na(domain_status) & !within_domain] <-
    "outside_expanded_simulation_domain"

n_eligible <- sum(eligible)
if (n_eligible < 2L) stop("Fewer than two eligible VMRs; score is undefined")
raw <- as.numeric(d[[score_basis]][eligible])
midrank <- rank(raw, ties.method = "average")
score <- (midrank - 0.5) / n_eligible
score_z <- as.numeric(scale(score))
if (any(!is.finite(score_z))) stop("Eligible score has zero or nonfinite variance")

d$local_snp_contribution_score_basis <- score_basis
d$local_genetic_control_eligible <- eligible
d$local_genetic_control_exclusion_reason <- reason
d$local_snp_contribution_score <- NA_real_
d$local_snp_contribution_score_z <- NA_real_
d$local_snp_contribution_quartile <- NA_character_
d$local_snp_contribution_score[eligible] <- score
d$local_snp_contribution_score_z[eligible] <- score_z
d$local_snp_contribution_quartile[eligible] <- ifelse(
    score <= 0.25, "bottom_quartile",
    ifelse(score >= 0.75, "top_quartile", "middle_50_percent")
)
d$absolute_pve_interpretation_allowed <- FALSE
d$local_genetic_control_decision <-
    "PASS_RELATIVE_GENETIC_CONTROL_FAIL_ABSOLUTE_LOCUS_PVE"
if ("positive_signal" %in% names(d)) {
    if ("positive_signal_audit_only" %in% names(d)) {
        stop("Input has both positive_signal and positive_signal_audit_only")
    }
    names(d)[names(d) == "positive_signal"] <- "positive_signal_audit_only"
}

write_tsv(d, cli$output)
cat("Wrote local SNP contribution scores for", n_eligible, "of", nrow(d),
    "VMRs in", unique(d$cohort), "x", unique(d$region), "\n")
