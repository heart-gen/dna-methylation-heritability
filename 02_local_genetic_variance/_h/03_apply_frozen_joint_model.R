#!/usr/bin/env Rscript

## Stage 03: apply the checksum-pinned final joint model once. The raw estimate
## is retained for audit; this stage does not authorize absolute-PVE language.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
h_dir <- dirname(script_path)
source(file.path(h_dir, "00_functions.R"))
source(file.path(h_dir, "joint_pve_functions.R"))

cli <- parse_cli(list(run_dir = ""))
if (!nzchar(cli$run_dir)) stop("--run-dir is required")
run_dir <- normalizePath(cli$run_dir)
manifest <- read_tsv(file.path(run_dir, "manifest.tsv"))
mval <- function(field) {
    value <- manifest$value[manifest$field == field]
    if (length(value) != 1L) stop("Run manifest lacks unique field: ", field)
    as.character(value[[1L]])
}
features <- read_tsv(file.path(
    run_dir, "results", "combined", "observed-joint-features.tsv"
))
model_path <- mval("joint_model_path")
development_path <- mval("development_features_path")
sha <- system2("sha256sum", normalizePath(model_path), stdout = TRUE)
observed_sha <- tolower(sub(" .*$", "", sha[[1L]]))
if (!identical(observed_sha, tolower(mval("joint_model_sha256")))) {
    stop("Frozen joint-model checksum changed after Stage 00")
}
model <- readRDS(model_path)
if (!identical(model$gate_version, "final_joint_pve_v1")) {
    stop("Unexpected joint-model gate version: ", model$gate_version)
}
development <- read_tsv(development_path)
development <- development[development$feature_complete %in% TRUE, , drop = FALSE]
domain_fields <- c("n", "num_snps", "p_eff", "ld_metric")
if (!all(domain_fields %in% names(development))) {
    stop("Development features lack domain fields")
}
domain <- list(
    num_snps = range(as.numeric(development$num_snps), finite = TRUE),
    ## p_eff was a measured feature, not a gridded design factor. Training
    ## minima/maxima are random order statistics, so the defensible support is
    ## its mathematical range rather than those accidental endpoints.
    ld_metric = range(as.numeric(development$ld_metric), finite = TRUE)
)
allowed_n <- sort(unique(as.integer(development$n)))

features$joint_pve_domain_status <- "feature_incomplete"
features$joint_pve_domain_reason <- ifelse(
    features$feature_complete %in% TRUE, NA_character_,
    "joint features incomplete"
)
complete <- which(features$feature_complete %in% TRUE &
                  !(features$computational_failure %in% TRUE))
if (length(complete)) {
    reasons <- vapply(complete, function(i) {
        reason <- character()
        if (!as.integer(features$n[[i]]) %in% allowed_n) {
            reason <- c(reason, "n outside simulated values")
        }
        p <- as.numeric(features$num_snps[[i]])
        if (!is.finite(p) || p < domain$num_snps[[1L]] ||
            p > domain$num_snps[[2L]]) {
            reason <- c(reason, "num_snps outside simulated range")
        }
        p_eff <- as.numeric(features$p_eff[[i]])
        if (!is.finite(p_eff) || p_eff < 1 ||
            p_eff > as.numeric(features$n[[i]])) {
            reason <- c(reason, "p_eff outside mathematical range [1,n]")
        }
        ld <- as.numeric(features$ld_metric[[i]])
        if (!is.finite(ld) || ld < domain$ld_metric[[1L]] ||
            ld > domain$ld_metric[[2L]]) {
            reason <- c(reason, "ld_metric outside development range")
        }
        if (length(reason)) paste(reason, collapse = "; ") else NA_character_
    }, character(1L))
    within <- is.na(reasons)
    features$joint_pve_domain_status[complete] <- ifelse(
        within, "within_domain", "outside_domain"
    )
    features$joint_pve_domain_reason[complete] <- reasons

    predicted <- predict_joint_pve(model, features[complete, , drop = FALSE])
    if (nrow(predicted) != length(complete)) stop("Prediction row-count mismatch")
    for (field in names(predicted)) features[[field]] <- NA
    for (field in names(predicted)) features[[field]][complete] <- predicted[[field]]
    if ("positive_signal" %in% names(features)) {
        names(features)[names(features) == "positive_signal"] <-
            "positive_signal_audit_only"
    }
} else {
    stop("No complete observed joint-feature rows")
}
features$joint_model_run_id <- mval("joint_model_run_id")
features$joint_model_sha256 <- observed_sha
features$absolute_pve_interpretation_allowed <- FALSE
write_tsv(features, file.path(
    run_dir, "results", "combined", "observed-joint-estimates.tsv"
))
cat("Applied frozen joint model to", length(complete), "VMRs\n")
