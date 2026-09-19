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
hyp <- cfg$axis$hypothesized_sign
main <- axis[outcome == cfg$axis$outcome]
row_for <- function(member) {
    ## A member is a spec (arm "base") or a primary-spec arm.
    r <- main[(spec == member & arm == "base") | (spec == "primary" & arm == member)]
    if (nrow(r) > 1L) stop("Ambiguous region-reading member: ", member)
    r
}
prim <- row_for("primary")
primary_supported <- prim$direction == hyp && prim$p < as.numeric(rr$primary_alpha)

member_rows <- list()
for (m in unlist(rr$same_n_members)) {
    r <- row_for(m)
    if (!nrow(r)) {
        ## Not fitted in this region (scMD fails its integration gate outside
        ## caudate). Recorded, and it does not count against the conjunction:
        ## an arm that cannot be fitted is absent evidence, not contrary evidence.
        member_rows[[m]] <- data.table(member = m, rule = rr$same_n_rule,
            fitted = FALSE, survives = NA, estimate = NA_real_, p = NA_real_)
        next
    }
    member_rows[[m]] <- data.table(member = m, rule = rr$same_n_rule, fitted = TRUE,
        survives = r$direction == hyp && r$p < as.numeric(rr$primary_alpha),
        estimate = r$estimate, p = r$p)
}
for (m in unlist(rr$reduced_n_members)) {
    r <- row_for(m)
    if (!nrow(r)) stop("Reduced-n member ", m, " was not fitted")
    frac <- r$estimate / prim$estimate
    member_rows[[m]] <- data.table(member = m, rule = rr$reduced_n_rule, fitted = TRUE,
        survives = r$direction == hyp && is.finite(frac) &&
            frac >= as.numeric(rr$reduced_n_min_fraction),
        estimate = r$estimate, p = r$p,
        fraction_of_primary = frac)
}
members <- rbindlist(member_rows, fill = TRUE)
all_survive <- all(members[fitted == TRUE, survives])

reading <- if (!primary_supported) {
    if (prim$p < as.numeric(rr$primary_alpha)) "OPPOSITE_DIRECTION" else "NOT_SUPPORTED"
} else if (all_survive) {
    "SUPPORTED_SURVIVES_GATING_SENSITIVITIES"
} else {
    "PRIMARY_ONLY_FAILS_GATING_SENSITIVITY"
}
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
