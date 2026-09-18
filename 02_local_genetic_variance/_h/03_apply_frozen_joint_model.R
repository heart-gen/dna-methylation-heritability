#!/usr/bin/env Rscript

## Stage 03: apply the checksum-pinned final joint model once. The raw estimate
## is retained for audit; this stage does not authorize absolute-PVE language.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
h_dir <- dirname(script_path)
source(file.path(h_dir, "00_functions.R"))
source(file.path(h_dir, "joint_pve_functions.R"))
## This stage verifies the support table against the sha Stage 00 pinned, so it
## needs file_sha256(). 00_shared/sha.R is dependency-free by design: the full
## loader pulls in config.R, which needs yaml, and yaml is deliberately absent
## from the calibrated-local-h2 env this stage runs in.
source(file.path(h_dir, "..", "..", "00_shared", "sha.R"))

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

## The domain is the support over which the frozen model has actually been
## characterised for THIS CELL: the union of the AR(1) training grid and only
## this cell's own observed-regime grids. Bounding p_eff only by its
## mathematical range [1, n] made this gate blind: it passed every eligible
## locus of lgv-AA-caudate-20260822 while 18.75% of them sat below the AR(1)
## minimum p_eff of 24.34, a regime the training grid never visited.
##
## Selecting the cell is not cosmetic. A support table pooled across cells lets
## one donor group widen another's eligibility domain -- admitting the EA cells
## dropped the shared p_eff floor from 2.058 to 1.344 -- and LD structure and
## MAF spectrum differ between donor groups, so that borrowing is not
## defensible. A table without a `cell` column is a pre-2026-09-13 pooled
## table; refuse it rather than silently reinstating the union.
support_path <- file.path(run_dir, "config",
                          "joint-pve-characterized-support.tsv")

## The support table is pinned the same way the model is. Stage 00 records
## config_characterized_support_sha256 from the repository copy and snapshots
## the file into the run; without this check the two could diverge -- an edited
## repository config, or an in-run copy replaced by hand -- and the run would
## still seal and still reach a README acceptance row citing a sha it no longer
## used. That is exactly the drift that happened once already, so it is a hard
## stop, not a warning. Stage 00 is the only place the field may be written.
support_sha <- tolower(file_sha256(support_path))
manifest_support_sha <- tolower(mval("config_characterized_support_sha256"))
if (is.na(support_sha)) {
    stop("Cannot checksum the characterized-support table: ", support_path)
}
if (!identical(support_sha, manifest_support_sha)) {
    stop("Characterized-support checksum changed after Stage 00.\n",
         "  manifest: ", manifest_support_sha, "\n",
         "  in-run:   ", support_sha, "\n",
         "  The eligibility domain is part of this run's provenance; refusing ",
         "to score against an unattributable support table (AGENTS.md 9).")
}

support <- read_tsv(support_path)
required_support <- c("cell", "feature", "support_min", "support_max",
                      "allowed_n")
if (!all(required_support %in% names(support))) {
    stop("Characterized-support table lacks required columns (needs a per-cell",
         " table with a `cell` column): ",
         paste(setdiff(required_support, names(support)), collapse = ", "))
}
run_cell <- mval("cohort")
support <- support[as.character(support$cell) == run_cell, , drop = FALSE]
if (!nrow(support)) {
    stop("Characterized-support table has no rows for cell: ", run_cell)
}
if (!all(c("num_snps", "p_eff", "ld_metric") %in% support$feature)) {
    stop("Characterized-support table lacks a required feature")
}
support_range <- function(field) {
    row <- support[support$feature == field, , drop = FALSE]
    if (nrow(row) != 1L) stop("Characterized support is not unique for ", field)
    values <- c(as.numeric(row$support_min[[1L]]),
                as.numeric(row$support_max[[1L]]))
    if (!all(is.finite(values))) stop("Nonfinite characterized support: ", field)
    values
}
domain <- list(
    num_snps = support_range("num_snps"),
    p_eff = support_range("p_eff"),
    ld_metric = support_range("ld_metric")
)
allowed_n <- sort(unique(as.integer(unlist(strsplit(
    as.character(support$allowed_n[[1L]]), ",", fixed = TRUE)))))
if (!length(allowed_n)) stop("Characterized support lacks allowed_n")

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
            reason <- c(reason, "n outside characterized values")
        }
        p <- as.numeric(features$num_snps[[i]])
        if (!is.finite(p) || p < domain$num_snps[[1L]] ||
            p > domain$num_snps[[2L]]) {
            reason <- c(reason, "num_snps outside characterized support")
        }
        p_eff <- as.numeric(features$p_eff[[i]])
        if (!is.finite(p_eff) || p_eff < 1 ||
            p_eff > as.numeric(features$n[[i]])) {
            reason <- c(reason, "p_eff outside mathematical range [1,n]")
        } else if (p_eff < domain$p_eff[[1L]] || p_eff > domain$p_eff[[2L]]) {
            reason <- c(reason, "p_eff outside characterized support")
        }
        ld <- as.numeric(features$ld_metric[[i]])
        if (!is.finite(ld) || ld < domain$ld_metric[[1L]] ||
            ld > domain$ld_metric[[2L]]) {
            reason <- c(reason, "ld_metric outside characterized support")
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
## Carried per row alongside the model sha, so the support that defined each
## locus's domain status is recoverable from the table alone.
features$config_characterized_support_sha256 <- support_sha
features$absolute_pve_interpretation_allowed <- FALSE
write_tsv(features, file.path(
    run_dir, "results", "combined", "observed-joint-estimates.tsv"
))
cat("Applied frozen joint model to", length(complete), "VMRs\n")
