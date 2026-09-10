#!/usr/bin/env Rscript
#### 10_environmental_exploratory -- coverage gate and interpretation flags ####
##
## Usage:
##   Rscript _h/04_apply_gates.R --run-id env-AA-caudate-YYYYMMDD
##
## This module's gate is deliberately NOT a scientific success criterion. A null
## axis result is a legitimate outcome and must not block sealing: AGENTS.md 2.3
## is explicit that "Low local SNP variance or poor SNP prediction is not
## evidence that a VMR is environmentally determined", and the converse holds
## here -- finding nothing along the axis says nothing about whether these VMRs
## respond to exposure. What the gate checks is whether the run had enough
## coverage to have said anything at all.
##
## The second job is to stamp the interpretation constraints onto a row, so a
## downstream consumer inherits them from the data rather than from prose in a
## README it may never open.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "10_environmental_exploratory", "_h")),
    "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "10_environmental_exploratory"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
env <- load_run_config("environmental", run_dir)

elig <- fread(file.path(run_dir, "results", "exposure-eligibility.tsv"))
assoc <- fread(file.path(run_dir, "results", "vmr-exposure-association.tsv"))
axis <- fread(file.path(run_dir, "results", "control-axis-test.tsv"))

n_eligible <- nrow(elig[eligible == TRUE])   # exposure x stratum pairs
n_tested_vmrs <- length(unique(assoc$vmr_id))
min_eligible <- as.integer(env$gates$min_eligible_exposures)
min_vmrs <- as.integer(env$gates$min_tested_vmrs)

checks <- data.table(
    check = c("eligible_exposures", "tested_vmrs", "axis_models_fitted",
              "no_forbidden_columns", "fdr_families_match_eligible"),
    observed = c(n_eligible, n_tested_vmrs, nrow(axis[status == "ok"]),
                 length(intersect(unlist(env$forbidden_columns),
                                  c(names(assoc), names(axis)))),
                 nrow(axis)),
    required = c(min_eligible, min_vmrs, 1L, 0L, n_eligible),
    comparison = c(">=", ">=", ">=", "==", "==")
)
checks[, pass := mapply(function(o, r, cmp) {
    switch(cmp, ">=" = o >= r, "==" = o == r, stop("bad comparison"))
}, observed, required, comparison)]

decision <- if (all(checks$pass)) "PASS_EXPLORATORY_COVERAGE" else
    "FAIL_EXPLORATORY_COVERAGE"

## A directional summary that cannot be mistaken for a claim: it names how many
## exposures moved along the axis, and immediately states what that does not
## license.
n_axis_fdr <- nrow(axis[status == "ok" & !is.na(primary_fdr) &
                        primary_fdr <= as.numeric(env$testing$fdr_alpha)])

dec <- data.table(
    run_id = opts$run_id, cohort = mf("cohort"), region = mf("region"),
    decision = decision,
    n_eligible_exposures = n_eligible,
    n_ineligible_exposures = nrow(elig[eligible == FALSE & kind == "binary"]),
    strata_tested = paste(sort(unique(axis$stratum)), collapse = ","),
    n_eligible_in_schizophrenia_stratum =
        nrow(elig[eligible == TRUE & stratum == "schizophrenia"]),
    n_near_gate_exposures = nrow(elig[near_gate == TRUE]),
    n_tested_vmrs = n_tested_vmrs,
    n_significant_vmr_exposure_pairs = sum(assoc$significant, na.rm = TRUE),
    n_axis_associations_fdr = n_axis_fdr,
    main_text_retention = "NEVER_SUPPLEMENT_ONLY",
    exploratory_supplement_only = TRUE,
    causal_interpretation_allowed = FALSE,
    absolute_pve_interpretation_allowed = FALSE,
    cross_region_comparison_allowed = FALSE,
    environmentally_determined_claim_allowed = FALSE,
    collider_flagged_exposures = paste(unlist(
        env$interpretation$collider_flagged_exposures), collapse = ","),
    burden_indicator_exposures = paste(unlist(
        env$interpretation$burden_indicator_exposures), collapse = ","),
    null_result_meaning = as.character(env$interpretation$null_result_meaning),
    cross_region_rationale = as.character(env$interpretation$cross_region_rationale)
)

write_atomic(checks, file.path(run_dir, "results", "gate-checks.tsv"))
write_atomic(dec, file.path(run_dir, "results", "environmental-decision.tsv"))

print(checks)
message("[10] decision: ", decision)
if (!all(checks$pass)) {
    stop("Coverage gate failed; see results/gate-checks.tsv. This is a ",
         "coverage failure, not a scientific finding.")
}
