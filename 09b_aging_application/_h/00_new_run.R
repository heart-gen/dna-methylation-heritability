#!/usr/bin/env Rscript
#### 09b_aging_application -- open a run ####
##
## Usage:
##   Rscript _h/00_new_run.R --cohort AA --region dlpfc [--allow-unlocked]
##
## Four upstreams: 01 for the per-VMR methylation phenotypes and donor
## covariates, 02 for the relative local-control score, 04 for the technical and
## chromatin covariates the axis model adjusts on, and 07 for the transcriptional
## coupling the descriptive annotation stage reads. All four must describe the
## same vmr_set_id.
##
## One diagnostic is run here, before anything is fitted, because it bounds what
## the module can say rather than being a result: VMRs were selected on residual
## SD after regressing out methylation PCs 1-5 that were computed WITHOUT age
## adjustment (01_vmr_catalog/_h/01_analyze.R). If a methPC tracks age, the
## catalog under-samples regions whose variability is mostly age-driven. The
## correlation is measured and written down; it is not corrected, since
## correcting would mean re-deriving the catalog.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "09b_aging_application"
MODULE_TAG <- "age"

opts <- parse_v2_args(require = c("cohort", "region"))
allow_unlocked <- isTRUE(opts$allow_unlocked)

cfg <- load_config("aging")
assert_locked(list(aging = cfg), allow_unlocked = allow_unlocked)

## AGENTS.md 3 and 7.2: the only admissible architecture predictor.
if (!identical(cfg$axis$predictor, "local_snp_contribution_score_z")) {
    stop("config/aging.yml axis.predictor must be local_snp_contribution_score_z. ",
         "Legacy predictability and absolute PVE are banned (AGENTS.md 3).")
}
## The headline's variance is the donor bootstrap plus the chromosome block
## jackknife (config header). Refuse a config that would thin the bootstrap
## below a usable variance estimate, or drop the diagnosis stratification.
n_boot_prod <- as.integer(cfg$inference$donor_bootstrap_n)
if (!is.finite(n_boot_prod) || n_boot_prod < 500L) {
    stop("inference.donor_bootstrap_n must be >= 500 for a production run.")
}
if (!identical(cfg$inference$bootstrap_stratified_by, "diagnosis")) {
    stop("inference.bootstrap_stratified_by must be 'diagnosis': cases are ~7 ",
         "years older, and every draw must keep the design the model adjusts for.")
}
if (!identical(cfg$axis$outcome, "debiased_sq_beta")) {
    stop("axis.outcome must be debiased_sq_beta. |beta|, its rank, |t| and ",
         "-log10 p are coupled to the score through the SE (config/aging.yml).")
}
for (flag in c("causal_interpretation_allowed",
               "environmentally_determined_claim_allowed",
               "absolute_pve_interpretation_allowed",
               "cross_region_raw_score_comparison_allowed")) {
    if (isTRUE(cfg$interpretation[[flag]])) {
        stop("interpretation.", flag, " must be false (AGENTS.md 2.3, 7.2, 8.1).")
    }
}
if (!identical(opts$cohort, cfg$cohort)) {
    stop("This module is prespecified for cohort ", cfg$cohort, ", not ", opts$cohort)
}
if (!opts$region %in% unlist(cfg$regions)) {
    stop("Region ", opts$region, " is not in config/aging.yml regions")
}

accepted <- lapply(cfg$upstreams, require_accepted_upstream,
                   cohort = opts$cohort, region = opts$region,
                   allow_unaccepted = allow_unlocked)

sets <- unlist(lapply(accepted, function(a) a$vmr_set_id))
sets <- unique(sets[!is.na(sets)])
if (length(sets) > 1) {
    stop("vmr_set_id mismatch across upstreams: ", paste(sets, collapse = " vs "),
         "\n  01, 02, 04 and 07 must describe the same VMR set (AGENTS.md 6).")
}
if (any(vapply(accepted, function(a) is.na(a$run_id %||% NA_character_), logical(1)))) {
    stop("An upstream has no accepted run for ", opts$cohort, " x ", opts$region,
         ". This module has no fallback to unaccepted runs: the stages read the ",
         "upstream run directories by ID.")
}

if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked); ",
         "production run IDs are derived, not chosen.")
}

n_boot <- if (allow_unlocked) as.integer(cfg$inference$smoke_donor_bootstrap_n) else
    n_boot_prod

run <- new_run(
    module = MODULE_TAG, cohort = opts$cohort, region = opts$region,
    module_root = file.path(repo_root(), MODULE),
    run_id = opts$run_id,
    vmr_set_id = if (length(sets)) sets[1] else NA_character_,
    upstream = stats::setNames(
        lapply(accepted, function(a) a$run_id %||% NA_character_),
        paste0(names(accepted), "_run_id")),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        config_aging_sha256 = file_sha256(file.path(repo_root(), "config", "aging.yml")),
        n_bootstrap = as.character(n_boot),
        axis_predictor = cfg$axis$predictor,
        axis_outcome = cfg$axis$outcome,
        axis_outcome_scale = cfg$axis$outcome_scale,
        bootstrap_stratified_by = cfg$inference$bootstrap_stratified_by
    )
)
for (d in c("results", "checkpoint", "logs")) {
    dir.create(file.path(run$dir, d), showWarnings = FALSE, recursive = TRUE)
}

## ------------------------------------------------ methPC-age catalog diagnostic
vmr_run_dir <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs",
                         accepted$vmr_catalog$run_id)
pc_files <- list.files(file.path(vmr_run_dir, "pca"), pattern = "^pc\\.csv$",
                       recursive = TRUE, full.names = TRUE)
if (!length(pc_files)) stop("No Module 01 methylation-PC tables under ", vmr_run_dir)
covar_prefix <- if (identical(opts$catalog_cohort, "AA")) "TOPMed_LIBD.AA" else "TOPMed_LIBD"
qcovar <- fread(file.path(vmr_run_dir, "covs", "chr_1", paste0(covar_prefix, ".qcovar")),
                header = FALSE, col.names = c("FID", "IID", "age"),
                colClasses = list(character = 1:2))

diag <- rbindlist(lapply(pc_files, function(f) {
    pc <- fread(f, colClasses = list(character = "brnum"))
    idx <- match(pc$brnum, qcovar$FID)
    if (anyNA(idx)) {
        stop("Methylation-PC donors missing from the Module 01 covariates: ", f)
    }
    chrom <- basename(dirname(f))
    pcs <- intersect(paste0("PC", 1:5), names(pc))
    rbindlist(lapply(pcs, function(p) data.table(
        chrom = chrom, pc = p,
        spearman_rho_age = stats::cor(pc[[p]], qcovar$age[idx], method = "spearman"),
        n = nrow(pc))))
}))
write_atomic(diag, file.path(run$dir, "results", "diagnostic-methpc-age.tsv"))
worst <- diag[which.max(abs(spearman_rho_age))]
append_manifest(run, list(
    methpc_age_max_abs_rho = sprintf("%.4f", abs(worst$spearman_rho_age)),
    methpc_age_max_abs_rho_at = paste0(worst$chrom, ":", worst$pc),
    methpc_age_median_abs_rho_pc1 = sprintf("%.4f",
        stats::median(abs(diag[pc == "PC1", spearman_rho_age])))
))
message("[09b] methPC-age diagnostic: max |rho| = ",
        sprintf("%.3f", abs(worst$spearman_rho_age)), " at ",
        worst$chrom, ":", worst$pc, " (written, not corrected)")

message("[09b] run ", run$run_id, " opened (", n_boot, " donor bootstrap draws)")
cat(run$run_id, "\n", sep = "")
