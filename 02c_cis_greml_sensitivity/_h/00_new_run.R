#!/usr/bin/env Rscript
#### 02c_cis_greml_sensitivity -- open a run ####
##
## Usage:
##   Rscript _h/00_new_run.R --cohort AA --region dlpfc [--allow-unlocked]
##
## The task universe is EXACTLY the accepted Module 02 run's per-VMR table --
## every VMR it scored or excluded -- so 02c fits the loci Module 02 saw and the
## comparison can never silently drop or add one. Module 02's score, its
## eligibility and its joint-feature estimators are copied into the task
## manifest read-only; nothing is written back to Module 02.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02c_cis_greml_sensitivity", "_h")),
                 "cgs_functions.R"))

opts <- parse_v2_args(require = c("cohort", "region"))
allow_unlocked <- isTRUE(opts$allow_unlocked)
cfg <- load_config("cis_greml_sensitivity")
assert_locked(list(cis_greml_sensitivity = cfg), allow_unlocked = allow_unlocked)
assert_cgs_interpretation(cfg)
gcta_v <- assert_gcta_pin(cfg)
thr <- load_config("thresholds")
if (!identical(opts$cohort, cfg$cohort)) stop("02c is prespecified for cohort ", cfg$cohort)
if (!opts$region %in% unlist(cfg$regions)) stop("region not in config regions")
if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked).")
}

up <- lapply(cfg$upstreams, require_accepted_upstream, cohort = opts$cohort,
             region = opts$region, allow_unaccepted = allow_unlocked)
if (any(vapply(up, function(u) is.na(u$run_id %||% NA_character_), logical(1)))) {
    stop("02c reads both upstream run directories by ID; each needs an accepted run.")
}
sets <- unique(vapply(up, function(u) u$vmr_set_id, character(1)))
if (length(sets) != 1L) stop("vmr_set_id differs between Modules 01 and 02: ",
                             paste(sets, collapse = " vs "))

lgc <- load_local_genetic_control(up$local_genetic_variance$run_id, opts$region,
                                  opts$cohort, eligible_only = FALSE)
if (!identical(unique(lgc$vmr_set_id), sets)) stop("Module 02 table vmr_set_id mismatch")
need <- c(cfg$summaries$score_column, "local_snp_contribution_score_z",
          unlist(cfg$summaries$module02_estimator_columns), "num_snps")
miss <- setdiff(need, names(lgc))
if (length(miss)) stop("Module 02 table lacks: ", paste(miss, collapse = ", "))

tasks <- lgc[, .(vmr_id, chrom, start, end,
                 module02_eligible = local_genetic_control_eligible %in% c(TRUE, "TRUE"),
                 module02_exclusion_reason = local_genetic_control_exclusion_reason,
                 module02_num_snps = num_snps,
                 local_snp_contribution_score, local_snp_contribution_score_z,
                 he_h2, bslmm_pve)]
tasks <- tasks[!toupper(sub("^chr", "", chrom)) %in% c("X", "Y")]
if (allow_unlocked) {
    el <- which(tasks$module02_eligible)
    k <- min(length(el), as.integer(cfg$tasks$smoke_n_vmrs))
    tasks <- tasks[el[unique(round(seq(1, length(el), length.out = k)))]]
}
setorder(tasks, chrom, start)
tasks[, task := ceiling(seq_len(.N) / as.integer(cfg$tasks$vmrs_per_task))]

modes <- rbindlist(lapply(cfg$reml_modes, function(m)
    data.table(reml_mode = m$id, role = m$role, args = paste(unlist(m$args), collapse = " "))))
if (sum(modes$role == "primary") != 1L) stop("Exactly one reml_mode must be primary")

design_n <- load_config("cohorts")$donor_counts[[opts$cohort]][[opts$region]]$design_n
run <- new_run(module = cfg$module_tag, cohort = opts$cohort, region = opts$region,
               module_root = file.path(repo_root(), MODULE), run_id = opts$run_id,
               vmr_set_id = sets,
               upstream = list(vmr_catalog_run_id = up$vmr_catalog$run_id,
                               local_genetic_variance_run_id = up$local_genetic_variance$run_id),
               extra = list(
                   smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
                   config_cis_greml_sensitivity_sha256 = file_sha256(
                       file.path(repo_root(), "config", "cis_greml_sensitivity.yml")),
                   gcta_binary = cfg$gcta$binary, gcta_version = gcta_v,
                   reml_primary_mode = modes[role == "primary", reml_mode],
                   design_n = as.character(design_n),
                   cis_window_bp = as.character(thr$cis$window_bp),
                   min_cis_variants = as.character(thr$cis$min_cis_variants),
                   n_vmrs = as.character(nrow(tasks)),
                   n_module02_eligible = as.character(sum(tasks$module02_eligible)),
                   replaces_module_02_score = "FALSE",
                   per_vmr_absolute_h2_reportable = "FALSE",
                   greml_significance_class_allowed = "FALSE"))
for (d in c("config", "results", "status", "summary", "figures", "work")) {
    dir.create(file.path(run$dir, d), recursive = TRUE, showWarnings = FALSE)
}
write_atomic(tasks, file.path(run$dir, "config", "vmr-tasks.tsv"))
write_atomic(modes, file.path(run$dir, "config", "reml-modes.tsv"))
append_manifest(run, list(n_tasks = max(tasks$task),
                          n_expected_units = nrow(tasks) * nrow(modes)))
cat(run$run_id, "\n", sep = "")
