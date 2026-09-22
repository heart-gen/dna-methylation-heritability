#!/usr/bin/env Rscript
#### 10_environmental_exploratory stage 06 -- collate the sealed runs ####
##
## Usage:
##   Rscript _h/06_collate_regions.R --cohort AA
##   Rscript _h/06_collate_regions.R --cohort AA --allow-unaccepted-runs
##   Rscript _h/06_collate_regions.R --cohort AA --allow-unlocked \
##       [--runs caudate=env-...,dlpfc=env-...,hippocampus=env-...]
##
## WHAT THIS IS NOT. It is not a cross-region stage in the sense of 09's stage 15
## or 09b's stage 05. `config/environmental.yml:interpretation` sets
## `cross_region_comparison_allowed: false`, and its `cross_region_rationale`
## says every result here is region-specific and no cross-region contrast is
## emitted: caudate is sequencing batch 3, region is perfectly confounded with
## batch (AGENTS.md 8.1), so a between-region exposure difference is
## uninterpretable. This stage therefore COLLATES and never COMPARES. It emits
## no pooled p, no region-general token, no between-region difference, and no
## decision of its own. It re-fits nothing.
##
## WHY IT EXISTS. The per-region results live in immutable `_m/runs/{RUN_ID}/`,
## which is gitignored and stays on Quest. A collaborator reading the manuscript
## needs the stage B table, the exposures that were testable, and the handful of
## stage A hits without re-running the SLURM chain. The tables written here are
## small, tracked in Git, and carry enough provenance to tie every number back
## to the run that produced it.
##
## FDR IS NOT RECOMPUTED. Every q-value is the one the sealed run wrote, from BH
## within that region's own family set. Pooling the three regions' p-values would
## be a cross-region operation and would silently change sealed numbers.
##
## ACCEPTANCE. Stages 17/18 of Module 09 stamp `-UNACCEPTED` on their filenames
## because their CONFIG is not PI-locked. That is not the situation here:
## `environmental.yml` is `pi_locked: true`. What is missing is the PI's
## acceptance row for the sealed runs. That is carried as data --
## `built_with_unaccepted_runs` and `citable` on every emitted table -- rather
## than as a filename, so accepting the runs and re-running this stage overwrites
## the same paths and the diff shows exactly the acceptance flip.
source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "10_environmental_exploratory"

## `parse_v2_args()` treats only `--allow-unlocked` as valueless, and it is
## shared by every module, so this valueless flag is stripped here rather than
## widening the shared parser for one stage's sake.
argv <- commandArgs(trailingOnly = TRUE)
allow_unaccepted_runs <- "--allow-unaccepted-runs" %in% argv
opts <- parse_v2_args(argv[argv != "--allow-unaccepted-runs"], require = "cohort")
allow_unlocked <- isTRUE(opts$allow_unlocked)
allow_unaccepted <- allow_unaccepted_runs || allow_unlocked

cfg <- load_config("environmental")
assert_locked(list(environmental = cfg), allow_unlocked = allow_unlocked)
interp <- cfg$interpretation
module_root <- file.path(repo_root(), MODULE)
out_dir <- file.path(module_root, "_m", "combined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

regions <- as.character(load_config("cohorts")$regions)
confounded <- as.character(config_get(load_config("region_donor_generalization"),
    "interpretation.technically_confounded_regions.all_outcomes"))
tier_of <- function(re) fifelse(re %in% confounded, "descriptive_only", "claim_eligible")

## ------------------------------------------------------------------- the runs
override <- list()
if (!is.null(opts$runs)) {
    if (!allow_unlocked) stop("--runs is for smoke use only (--allow-unlocked)")
    for (kv in strsplit(opts$runs, ",", fixed = TRUE)[[1]]) {
        p <- strsplit(kv, "=", fixed = TRUE)[[1]]
        override[[p[1]]] <- p[2]
    }
}
resolve_run <- function(region) {
    if (!is.null(override[[region]])) return(override[[region]])
    acc <- tryCatch(require_accepted_upstream(MODULE, opts$cohort, region)$run_id,
                    error = function(e) NA_character_)
    if (!is.na(acc)) return(acc)
    if (!allow_unaccepted) {
        stop("No accepted ", MODULE, " run for ", opts$cohort, " x ", region,
             ". Record it in ", MODULE, "/README.md under 'Accepted runs' ",
             "(AGENTS.md 6), or pass --allow-unaccepted-runs to collate the ",
             "sealed runs with citable = FALSE.", call. = FALSE)
    }
    cand <- list.files(file.path(module_root, "_m", "runs"),
                       pattern = paste0("^env-", opts$cohort, "-", region, "-"))
    if (!length(cand)) stop("No stage 10 run found for ", region)
    sort(cand, decreasing = TRUE)[[1L]]
}
runs <- vapply(regions, resolve_run, character(1))
run_dir <- function(re) file.path(module_root, "_m", "runs", runs[[re]])
res <- function(re, f) file.path(run_dir(re), "results", f)

manifests <- lapply(regions, function(re) {
    m <- fread(file.path(run_dir(re), "manifest.tsv"), colClasses = "character")
    setNames(as.list(m$value), m$field)
})
names(manifests) <- regions
mf <- function(re, f) {
    v <- manifests[[re]][[f]]
    if (is.null(v)) NA_character_ else v
}

## A sealed run is the precondition for collating it at all: an open run can
## still change under us, and then a tracked table would describe a state that
## no longer exists.
unsealed <- regions[vapply(regions, function(re) !nzchar(mf(re, "sealed_at")) ||
                               is.na(mf(re, "sealed_at")), logical(1))]
if (length(unsealed)) {
    stop("Not sealed, so not collatable: ", paste(runs[unsealed], collapse = ", "),
         ". Run _h/05_finalize_run.R first.", call. = FALSE)
}
smoke <- any(vapply(regions, function(re) identical(mf(re, "smoke_run"), "TRUE"),
                    logical(1)))

## Upstream currency, region by region. A run built on a Module 02 score that is
## no longer the accepted one is still readable, but a writer needs to know.
upstream <- rbindlist(lapply(regions, function(re) {
    cited <- mf(re, "upstream_local_genetic_variance_run_id")
    acc <- tryCatch(
        require_accepted_upstream("02_local_genetic_variance", opts$cohort, re)$run_id,
        error = function(e) NA_character_)
    data.table(region = re, run_id = runs[[re]],
               module_02_cited = cited, module_02_accepted = acc,
               upstream_current = identical(cited, acc))
}))
n_stale <- sum(!upstream$upstream_current)

citable <- !allow_unaccepted && !smoke && n_stale == 0L
prov_cols <- function(dt, re) {
    dt[, `:=`(cohort = opts$cohort, region = re, run_id = runs[[re]],
              tier = tier_of(re), sealed_at = mf(re, "sealed_at"),
              run_git_commit = mf(re, "git_commit"),
              config_environmental_sha256 = mf(re, "config_environmental_sha256"),
              vmr_set_id = mf(re, "vmr_set_id"),
              citable = citable, built_with_unaccepted_runs = allow_unaccepted,
              cross_region_comparison_allowed = FALSE,
              cross_region_contrast_emitted = FALSE,
              regions_are_independent_replicates = FALSE,
              exploratory_supplement_only = TRUE)]
    dt
}
stack <- function(f, ...) rbindlist(lapply(regions, function(re) {
    p <- res(re, f)
    if (!file.exists(p)) return(NULL)
    dt <- fread(p, ...)
    ## The per-region tables already carry cohort/region/run_id on most rows;
    ## setting them again keeps every emitted table uniform and self-describing.
    prov_cols(dt, re)
}), fill = TRUE)

## ------------------------------------------------------- stage B, the headline
axis <- stack("control-axis-test.tsv")
setorder(axis, region, stratum, primary_fdr, na.last = TRUE)

## The 81-column run table is not readable by hand. This is a strict COLUMN
## SUBSET of it -- no new numbers, nothing recomputed -- holding what a writer
## needs to state a result and its two standing caveats.
read_cols <- c("region", "tier", "exposure", "stratum", "n_vmrs_in_axis",
               "n_donors_refit", "primary_scale", "primary_beta", "primary_se",
               "primary_ci_lower", "primary_ci_upper", "primary_p", "primary_fdr",
               "absolute_beta", "absolute_p", "absolute_role",
               "mean_omega", "mean_omega_z", "fieller_bounded", "fieller_role",
               "bootstrap_mean_omega_inflation", "arm_covariates", "arm_beta",
               "arm_p", "arm_attenuation", "neglog10p_beta", "neglog10p_p",
               "mechanically_biased_toward_hypothesis", "run_id", "citable",
               "built_with_unaccepted_runs")
axis_reading <- axis[, intersect(read_cols, names(axis)), with = FALSE]

## ------------------------------------------------------- stage A and the gates
families <- stack("fdr-families.tsv")
eligibility <- stack("exposure-eligibility.tsv")
gates <- stack("gate-checks.tsv")
decisions <- stack("environmental-decision.tsv")
decisions <- merge(decisions, upstream[, .(region, module_02_cited,
                                           module_02_accepted, upstream_current)],
                   by = "region", all.x = TRUE)

## Only the FDR-surviving per-VMR rows. The full per-VMR table is 11 MB a region
## and stays on Quest; these are the rows a manuscript can name.
hits <- rbindlist(lapply(regions, function(re) {
    dt <- fread(res(re, "vmr-exposure-association.tsv"))
    prov_cols(dt[significant == TRUE], re)
}), fill = TRUE)
if (nrow(hits)) setorder(hits, region, exposure, stratum, fdr)

## ------------------------------------------------------------------ provenance
provenance <- rbindlist(lapply(regions, function(re) data.table(
    cohort = opts$cohort, region = re, run_id = runs[[re]],
    tier = tier_of(re),
    sealed_at = mf(re, "sealed_at"),
    decision = mf(re, "decision"),
    vmr_set_id = mf(re, "vmr_set_id"),
    n_donors = as.integer(mf(re, "n_donors_exposure_matrix")),
    donor_checksum = mf(re, "donor_checksum_exposure_matrix"),
    n_tested_vmrs = as.integer(mf(re, "n_tested_vmrs")),
    n_fdr_families = as.integer(mf(re, "n_fdr_families")),
    n_significant_pairs = as.integer(mf(re, "n_significant_pairs")),
    n_axis_associations_fdr = as.integer(mf(re, "n_axis_associations_fdr")),
    n_absolute_p_below_alpha = as.integer(mf(re, "n_absolute_p_below_alpha")),
    n_fieller_bounded = as.integer(mf(re, "n_fieller_bounded")),
    max_bootstrap_mean_omega_inflation =
        as.numeric(mf(re, "max_bootstrap_mean_omega_inflation")),
    axis_primary_estimand = mf(re, "primary_axis_estimand"),
    upstream_local_genetic_variance_run_id =
        mf(re, "upstream_local_genetic_variance_run_id"),
    upstream_repeat_architecture_run_id =
        mf(re, "upstream_repeat_architecture_run_id"),
    upstream_vmr_catalog_run_id = mf(re, "upstream_vmr_catalog_run_id"),
    run_git_commit = mf(re, "git_commit"),
    config_environmental_sha256 = mf(re, "config_environmental_sha256"))))
provenance <- merge(provenance,
                    upstream[, .(region, module_02_accepted, upstream_current)],
                    by = "region", all.x = TRUE)
provenance[, `:=`(
    collated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    collated_at_git_commit = git_commit(),
    collated_at_git_dirty = git_dirty(),
    citable = citable,
    built_with_unaccepted_runs = allow_unaccepted,
    smoke_run = smoke,
    n_regions_collated = length(regions),
    fdr_recomputed_across_regions = FALSE,
    cross_region_comparison_allowed = FALSE,
    cross_region_contrast_emitted = FALSE,
    regions_are_independent_replicates = FALSE,
    pooled_p_emitted = FALSE,
    exploratory_supplement_only = TRUE,
    causal_interpretation_allowed = FALSE,
    absolute_pve_interpretation_allowed = FALSE,
    environmentally_determined_claim_allowed = FALSE,
    null_result_meaning = as.character(interp$null_result_meaning),
    variance_budget_limitation = as.character(interp$variance_budget_limitation),
    cross_region_rationale = as.character(interp$cross_region_rationale))]

## ----------------------------------------------------------------------- write
sfx <- if (allow_unlocked) paste0("-", opts$cohort, "-UNLOCKED") else
    paste0("-", opts$cohort)
w <- function(x, stem) {
    if (is.null(x) || !nrow(x)) return(invisible(NULL))
    write_atomic(x, file.path(out_dir, paste0(stem, sfx, ".tsv")))
}
w(axis, "environmental-axis-per-region")
w(axis_reading, "environmental-axis-reading")
w(families, "environmental-fdr-families")
w(eligibility, "environmental-exposure-eligibility")
w(gates, "environmental-gate-checks")
w(decisions, "environmental-decision")
w(hits, "environmental-vmr-associations-fdr")
w(provenance, "environmental-collation-provenance")

cat(sprintf("Collated %d region(s) into %s\n", length(regions), out_dir))
for (re in regions) {
    cat(sprintf("  %-12s %s  sealed %s  tier %s  upstream_current %s\n",
                re, runs[[re]], mf(re, "sealed_at"), tier_of(re),
                upstream[region == re]$upstream_current))
}
cat(sprintf("  stage B families %d | stage A FDR-surviving VMR-exposure pairs %d\n",
            nrow(axis), nrow(hits)))
cat(sprintf("  citable = %s%s\n", citable,
            if (!citable) paste0("  (",
                paste(c(if (allow_unaccepted) "runs not accepted",
                        if (smoke) "smoke run",
                        if (n_stale) paste0(n_stale, " stale upstream")),
                      collapse = "; "), ")") else ""))
cat("  no cross-region contrast emitted (cross_region_comparison_allowed = FALSE)\n")
