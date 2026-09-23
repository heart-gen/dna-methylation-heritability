#!/usr/bin/env Rscript
#### 09b_aging_application -- coverage gate, region reading, interpretation flags ####
##
## Usage:
##   Rscript _h/03_apply_gates.R --run-id age-AA-dlpfc-YYYYMMDD
##
## Two separate things, kept separate:
##
##   1. The COVERAGE gate decides whether the run may seal. It is not a success
##      criterion; a null axis result is a legitimate outcome and seals.
##   2. The REGION READING states what the run supports, under strict
##      conjunction over the gating specs and arms (config/aging.yml:
##      region_reading). It never blocks sealing.
##
## Interpretation constraints are stamped onto the decision row so a downstream
## consumer inherits them from the data, not from a README it may not open.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
H_DIR <- Sys.getenv("V2_RUN_CODE",
                    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09b_aging_application", "_h"))
source(file.path(H_DIR, "run_config.R"))
## region_reading_members() lives there so a test can drive the reading from a
## sealed run's axis-tests.tsv without rerunning a stage.
source(file.path(H_DIR, "age_functions.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "09b_aging_application"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)
manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
if (run_is_sealed(manifest)) stop("Run is sealed and immutable: ", opts$run_id)
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mf("cohort"); region <- mf("region")
cfg <- load_run_config("aging", run_dir)
rr <- cfg$region_reading

axis <- fread(file.path(run_dir, "results", "axis-tests.tsv"))
eff <- fread(file.path(run_dir, "results", "vmr-age-effects.tsv"))
spec_summary <- fread(file.path(run_dir, "results", "age-spec-summary.tsv"))

## ------------------------------------------------------------- coverage gate
n_universe <- as.integer(mf("n_vmrs_universe"))
n_modelled <- as.integer(mf("n_vmrs_modelled"))
B <- as.integer(mf("n_bootstrap"))
checks <- data.table(
    check = c("modelled_vmrs", "nonfinite_vmr_fraction", "bootstrap_completed",
              "bootstrap_failed_fraction",
              "primary_spec_fitted", "no_forbidden_columns"),
    observed = c(n_modelled,
                 (n_universe - n_modelled) / n_universe,
                 as.integer(mf("n_bootstrap_completed")),
                 max(axis$n_bootstrap_failed) / B,
                 as.integer(nrow(axis[spec == "primary" & arm == "base" &
                                      outcome == cfg$axis$outcome]) == 1L),
                 length(intersect(unlist(cfg$forbidden_columns),
                                  c(names(axis), names(eff))))),
    required = c(as.numeric(cfg$gates$min_modelled_vmrs),
                 as.numeric(cfg$gates$max_nonfinite_vmr_fraction),
                 B, as.numeric(cfg$gates$max_failed_bootstrap_fraction), 1, 0),
    comparison = c(">=", "<=", "==", "<=", "==", "==")
)
checks[, pass := mapply(function(o, r, cmp) {
    switch(cmp, ">=" = o >= r, "<=" = o <= r, "==" = o == r, stop("bad comparison"))
}, observed, required, comparison)]
decision <- if (all(checks$pass)) "PASS_AGING_AXIS_COVERAGE" else "FAIL_AGING_AXIS_COVERAGE"

## ------------------------------------------------------------- region reading
## The conjunction is over FITTED members only, and a member can fail to be
## fitted in two ways, both of them the scMD integration gate: the `cell_scmd`
## SPEC is declined in 01_age_effects.R, and an scMD-derived gating ARM is
## declined in 02_axis_test.R. Both record a reason, collected here so the
## decision row says why rather than only that.
hyp <- cfg$axis$hypothesized_sign
main <- axis[outcome == cfg$axis$outcome]
skipped_f <- file.path(run_dir, "results", "axis-arms-skipped.tsv")
skipped_arms <- if (file.exists(skipped_f)) fread(skipped_f) else
    data.table(arm = character(0), reason = character(0))
pick <- function(tab, key_col, key) {
    if (!all(c(key_col, "reason") %in% names(tab))) return(character(0))
    v <- as.character(tab[["reason"]][tab[[key_col]] == key])
    v[!is.na(v)]
}
reason_for <- function(m) {
    r <- c(pick(spec_summary, "spec", m), pick(skipped_arms, "arm", m))
    if (length(r)) r[1] else NA_character_
}
rd <- region_reading_members(main, rr, hyp, reason_for)
members <- rd$members
prim <- rd$prim
primary_supported <- rd$primary_supported
all_survive <- rd$all_survive
reading <- rd$reading
region_supported <- reading == "SUPPORTED_SURVIVES_GATING_SENSITIVITIES"

confounded <- as.character(config_get(load_config("region_donor_generalization"),
    "interpretation.technically_confounded_regions.all_outcomes"))

dec <- data.table(
    run_id = opts$run_id, cohort = cohort, region = region,
    decision = decision,
    region_reading = reading,
    region_supported = region_supported,
    primary_estimate = prim$estimate,
    primary_estimate_meaning = paste0(
        "proportional change in mean squared age effect (debiased) per SD of ",
        "local_snp_contribution_score_z"),
    primary_se = prim$se,
    primary_ci_lower = prim$ci_lower,
    primary_ci_upper = prim$ci_upper,
    primary_p = prim$p,
    primary_direction = prim$direction,
    hypothesized_direction = hyp,
    n_vmrs_modelled = n_modelled,
    n_donors = as.integer(mf("n_donors")),
    n_bootstrap = B,
    gating_members_fitted = paste(members[fitted == TRUE, member], collapse = ","),
    gating_members_failing = paste(members[fitted == TRUE & survives == FALSE, member],
                                   collapse = ","),
    gating_members_not_fitted = paste(members[fitted == FALSE, member], collapse = ","),
    gating_members_not_fitted_reason = paste(
        members[fitted == FALSE, paste0(member, ":", reason)], collapse = ","),
    ## Which regions the scMD-derived members could be fitted in at all.
    scmd_integration_gate = mf("scmd_integration_gate"),
    technically_confounded_region = region %in% confounded,
    cross_sectional_design = TRUE,
    causal_interpretation_allowed = FALSE,
    environmentally_determined_claim_allowed = FALSE,
    absolute_pve_interpretation_allowed = FALSE,
    cross_region_raw_score_comparison_allowed = FALSE,
    epigenetic_clock_claim_allowed = FALSE,
    wording = as.character(cfg$interpretation$wording),
    catalog_scope_note = as.character(cfg$interpretation$catalog_scope_note),
    null_result_meaning = as.character(cfg$interpretation$null_result_meaning),
    methpc_age_max_abs_rho = as.numeric(mf("methpc_age_max_abs_rho")),
    smoke_run = identical(mf("smoke_run"), "TRUE")
)

write_atomic(checks, file.path(run_dir, "results", "gate-checks.tsv"))
write_atomic(members, file.path(run_dir, "results", "gating-sensitivities.tsv"))
write_atomic(dec, file.path(run_dir, "results", "aging-decision.tsv"))

print(checks)
print(members)
message("[09b] ", decision, "; region reading: ", reading)
if (!all(checks$pass)) {
    stop("Coverage gate failed; see results/gate-checks.tsv. This is a ",
         "coverage failure, not a scientific finding.")
}
