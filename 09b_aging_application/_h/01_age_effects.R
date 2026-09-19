#!/usr/bin/env Rscript
#### 09b_aging_application -- per-VMR age effects ####
##
## Usage:
##   Rscript _h/01_age_effects.R --run-id age-AA-dlpfc-YYYYMMDD
##
## Fits `meth ~ age + sex + diagnosis` (and the sensitivity specs of
## config/aging.yml:age_model.specs) for every Module 02-eligible VMR. The locus
## universe is Module 02's eligible loci, so every fitted VMR has a score, and
## the covariates come from the reader Module 02 used, so the age effect and the
## score condition on the same things.
##
## Everything is fitted as one matrix problem per spec, and the donor x VMR
## matrix is checkpointed: stages 02 and 05 bootstrap donors against it a
## thousand times.
##
## The per-VMR BH q is descriptive only. The module's inference is the axis
## test in stage 02, not a count of age-significant VMRs.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
H_DIR <- Sys.getenv("V2_RUN_CODE",
                    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09b_aging_application", "_h"))
source(file.path(H_DIR, "run_config.R"))
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

## ------------------------------------------------------------- locus universe
score <- load_local_genetic_control(mf("upstream_local_genetic_variance_run_id"),
                                    region = region, cohort = cohort)
banned <- intersect(unlist(cfg$forbidden_columns), names(score))
if (length(banned)) {
    stop("Module 02 table carries superseded column(s): ",
         paste(banned, collapse = ", "), " (AGENTS.md 3)")
}
audit <- intersect(unlist(cfg$audit_only_columns), names(score))
if (length(audit)) score[, (audit) := NULL]
if (!identical(unique(score$vmr_set_id), mf("vmr_set_id"))) {
    stop("Module 02 vmr_set_id does not match the run's")
}
tasks <- score[order(chrom, start), .(vmr_id, chrom, start, end, n_cpgs)]
lgv_n <- unique(as.integer(score$n))
if (length(lgv_n) != 1L) {
    stop("Module 02 eligible loci were fitted on differing donor counts: ",
         paste(lgv_n, collapse = ","))
}

## ------------------------------------------------------------- phenotypes
vmr_run_dir <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs",
                         mf("upstream_vmr_catalog_run_id"))
covar_prefix <- if (identical(cohort, "AA")) "TOPMed_LIBD.AA" else "TOPMed_LIBD"
message("[09b] reading ", nrow(tasks), " VMR phenotypes")
inp <- load_age_inputs(tasks, vmr_run_dir, catalog_cohort = cohort,
                       covar_prefix = covar_prefix)
Y <- inp$Y; donors <- inp$donors
if (!identical(colnames(Y), tasks$vmr_id)) stop("VMR column order lost")

## The design n must match what Module 02 fitted and what cohorts.yml locks.
if (nrow(donors) != lgv_n) {
    stop("Donor count ", nrow(donors), " differs from the Module 02 fits (",
         lgv_n, ") for the same VMRs")
}
assert_expected_n(nrow(donors), cohort, region)

## A VMR with any non-finite donor value cannot enter a shared-design matrix
## fit. It is excluded from every spec and counted; the gate caps the fraction.
finite_vmr <- colSums(!is.finite(Y)) == 0L
excluded <- tasks[!finite_vmr, .(vmr_id, reason = "nonfinite_phenotype")]
Y <- Y[, finite_vmr, drop = FALSE]
tasks <- tasks[finite_vmr]

## ------------------------------------------------------------- cell PCs
cc <- cfg$cell_composition
cellpcs <- list(
    music = cell_composition_pcs(
        file.path(repo_root(), sub("{region}", region, cc$rna_music_template, fixed = TRUE)),
        n_pcs = cc$n_pcs))
scmd_ok <- scmd_gate_passes(region, cc$scmd_gate_table)
if (scmd_ok) {
    cellpcs$scmd <- cell_composition_pcs(
        file.path(repo_root(), sub("{region}", region, cc$dnam_scmd_template, fixed = TRUE)),
        n_pcs = cc$n_pcs)
}

## ------------------------------------------------------------- fit every spec
specs <- cfg$age_model$specs
designs <- list()
spec_rows <- list()
effects <- list()
for (nm in names(specs)) {
    sp <- specs[[nm]]
    if (isTRUE(sp$requires_scmd_integration_gate) && !scmd_ok) {
        spec_rows[[nm]] <- data.table(spec = nm, role = sp$role, fitted = FALSE,
            reason = "scmd_integration_gate_fails_in_region",
            n_donors = NA_integer_, n_controls = NA_integer_, n_cases = NA_integer_)
        next
    }
    des <- build_spec_design(sp, donors, cfg, cellpcs)
    fit <- fit_age_matrix(des$X, Y[des$rows, , drop = FALSE], des$age_col)
    designs[[nm]] <- des
    spec_rows[[nm]] <- data.table(spec = nm, role = sp$role, fitted = TRUE,
        reason = NA_character_, n_donors = des$n,
        n_controls = des$n_controls, n_cases = des$n_cases,
        design_terms = paste(colnames(des$X), collapse = ","), df = fit$df)
    effects[[nm]] <- data.table(
        vmr_id = tasks$vmr_id, spec = nm, n_donors = des$n,
        beta_age_per_decade = fit$beta, se = fit$se, t = fit$t, p = fit$p,
        q_bh_descriptive = stats::p.adjust(fit$p, method = "BH"))
    message(sprintf("[09b] %-14s n=%d  median |beta|=%.4g per decade",
                    nm, des$n, stats::median(abs(fit$beta))))
}
if (is.null(designs$primary)) stop("The primary spec was not fitted")
if (designs$primary$n != nrow(donors)) {
    stop("The primary spec must use every donor; it used ", designs$primary$n)
}

spec_summary <- rbindlist(spec_rows, fill = TRUE)
eff <- rbindlist(effects)
eff[, `:=`(cohort = cohort, region = region, run_id = opts$run_id,
           cross_sectional_design = TRUE)]

write_atomic(eff, file.path(run_dir, "results", "vmr-age-effects.tsv"))
write_atomic(spec_summary, file.path(run_dir, "results", "age-spec-summary.tsv"))
write_atomic(excluded, file.path(run_dir, "results", "excluded-vmrs.tsv"))
saveRDS(list(Y = Y, donors = donors, vmr_ids = tasks$vmr_id,
             chrom = as.character(tasks$chrom),
             designs = designs, region = region, run_id = opts$run_id),
        file.path(run_dir, "checkpoint", "age-inputs.rds"))

append_manifest(list(dir = run_dir), list(
    n_donors = as.character(nrow(donors)),
    donor_checksum = donor_checksum(donors$FID),
    n_vmrs_universe = as.character(length(finite_vmr)),
    n_vmrs_modelled = as.character(ncol(Y)),
    n_vmrs_nonfinite = as.character(sum(!finite_vmr)),
    scmd_integration_gate = if (scmd_ok) "PASS" else "FAIL",
    specs_fitted = paste(names(designs), collapse = ",")
))
print(spec_summary)
