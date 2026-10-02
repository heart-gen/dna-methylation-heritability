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

## Upstream currency, region by region, over EVERY upstream the run pinned. A run
## built on an upstream that is no longer the accepted one is still readable, but
## a writer needs to know, and `citable` below turns on this.
##
## WIDENED 2026-10-02. This checked Module 02 alone, and `00_new_run.R` pins
## three. The gap was not hypothetical: env-AA-*-20260920-a pinned
## rra-AA-*-20260906, Module 04 superseded it with rra-AA-*-20260925-a on
## 2026-09-25, and `upstream_current` went on reporting TRUE on all three rows
## while `citable` stayed TRUE -- which is exactly the state this column exists
## to make visible. `require_accepted_upstream()` does no transitive check, so
## Module 11 would have consumed it silently too.
UPSTREAM_MODULES <- c(
    vmr_catalog            = "01_vmr_catalog",
    local_genetic_variance = "02_local_genetic_variance",
    repeat_architecture    = "04_repeat_repressive_architecture"
)
upstream <- rbindlist(lapply(regions, function(re) {
    per <- lapply(names(UPSTREAM_MODULES), function(key) {
        cited <- mf(re, paste0("upstream_", key, "_run_id"))
        acc <- tryCatch(
            require_accepted_upstream(UPSTREAM_MODULES[[key]], opts$cohort,
                                      re)$run_id,
            error = function(e) NA_character_)
        ## A run that pinned nothing for an upstream cannot be judged stale on
        ## it; that is a provenance hole, reported as NA rather than as current.
        list(key = key, cited = cited, accepted = acc,
             current = if (!nzchar(cited %||% "") || is.na(cited)) NA
                       else identical(cited, acc))
    })
    stale <- vapply(per, function(x) identical(x$current, FALSE), logical(1))
    unknown <- vapply(per, function(x) is.na(x$current), logical(1))
    row <- data.table(
        region = re, run_id = runs[[re]],
        ## Kept under their original names: these three columns are read by
        ## `environmental-decision-{cohort}.tsv` consumers written before today.
        module_02_cited = mf(re, "upstream_local_genetic_variance_run_id"),
        module_02_accepted = per[[which(names(UPSTREAM_MODULES) ==
                                        "local_genetic_variance")]]$accepted,
        n_upstreams_checked = length(per),
        n_upstreams_stale = sum(stale),
        n_upstreams_unpinned = sum(unknown),
        stale_upstreams = if (any(stale)) paste(vapply(per[stale], function(x)
            sprintf("%s: cited %s, accepted %s", x$key, x$cited, x$accepted),
            character(1)), collapse = "; ") else NA_character_,
        upstream_current = !any(stale) && !any(unknown))
    for (x in per) {
        data.table::set(row, j = paste0("upstream_", x$key, "_cited"),
                        value = x$cited)
        data.table::set(row, j = paste0("upstream_", x$key, "_accepted"),
                        value = x$accepted)
        data.table::set(row, j = paste0("upstream_", x$key, "_current"),
                        value = x$current)
    }
    row
}), fill = TRUE)
n_stale <- sum(!upstream$upstream_current)
if (n_stale) {
    for (i in which(!upstream$upstream_current)) {
        warning("Run ", upstream$run_id[i], " (", upstream$region[i],
                ") is not current: ",
                upstream$stale_upstreams[i] %||% "an upstream is unpinned",
                call. = FALSE)
    }
}

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

## The non-gating arms in long form, one row per region x exposure x stratum x
## arm. Written by stage 03 only when at least one arm is declared, so a run
## sealed before 2026-10-02 simply contributes nothing here.
arms_long <- stack("control-axis-arms.tsv")
if (nrow(arms_long)) setorder(arms_long, region, arm, stratum, fdr, na.last = TRUE)

## The 81-column run table is not readable by hand. This is a strict COLUMN
## SUBSET of it -- no new numbers, nothing recomputed -- holding what a writer
## needs to state a result and its two standing caveats.
read_cols <- c("region", "tier", "exposure", "stratum", "n_vmrs_in_axis",
               "n_donors_refit", "primary_scale", "primary_beta", "primary_se",
               "primary_ci_lower", "primary_ci_upper", "primary_p", "primary_fdr",
               "absolute_beta", "absolute_p", "absolute_role",
               "mean_omega", "mean_omega_z", "fieller_bounded", "fieller_role",
               "bootstrap_mean_omega_inflation", "n_axis_arms",
               "axis_arm_names", "neglog10p_beta", "neglog10p_p",
               "mechanically_biased_toward_hypothesis", "run_id", "citable",
               "built_with_unaccepted_runs")
## Arm columns are selected by PATTERN, not by name, so declaring a new arm in
## config does not silently drop it out of the reading table -- which is what a
## fixed list of `arm_beta`/`arm_p`/`arm_attenuation` did when the arms became a
## map on 2026-10-02.
arm_read <- grep("^arm_.*_(beta|p|fdr|attenuation|n_vmrs_lost_vs_primary)$",
                 names(axis), value = TRUE)
axis_reading <- axis[, c(intersect(read_cols, names(axis)), arm_read),
                     with = FALSE]

## ------------------------------------------------------- stage A and the gates
families <- stack("fdr-families.tsv")
eligibility <- stack("exposure-eligibility.tsv")
gates <- stack("gate-checks.tsv")
decisions <- stack("environmental-decision.tsv")
decisions <- merge(decisions, upstream[, setdiff(names(upstream), "run_id"),
                                       with = FALSE],
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
                    upstream[, c("region", "module_02_accepted",
                                 "upstream_current", "n_upstreams_checked",
                                 "n_upstreams_stale", "n_upstreams_unpinned",
                                 "stale_upstreams"), with = FALSE],
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
w(arms_long, "environmental-axis-arms")
w(families, "environmental-fdr-families")
w(eligibility, "environmental-exposure-eligibility")
w(gates, "environmental-gate-checks")
w(decisions, "environmental-decision")
w(hits, "environmental-vmr-associations-fdr")
w(provenance, "environmental-collation-provenance")

cat(sprintf("Collated %d region(s) into %s\n", length(regions), out_dir))
for (re in regions) {
    u <- upstream[region == re]
    cat(sprintf("  %-12s %s  sealed %s  tier %s  upstream_current %s (%d/%d)\n",
                re, runs[[re]], mf(re, "sealed_at"), tier_of(re),
                u$upstream_current,
                u$n_upstreams_checked - u$n_upstreams_stale -
                    u$n_upstreams_unpinned,
                u$n_upstreams_checked))
    if (!isTRUE(u$upstream_current)) {
        cat("                 stale: ",
            u$stale_upstreams %||% "an upstream is unpinned", "\n", sep = "")
    }
}
if (nrow(arms_long)) {
    cat(sprintf("  non-gating arms: %s\n",
                paste(sort(unique(arms_long$arm)), collapse = ", ")))
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
