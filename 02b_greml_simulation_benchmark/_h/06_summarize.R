#!/usr/bin/env Rscript
#### 02b -- reconcile every unit, then compute recovery metrics ####
##
## Reconciliation comes first and is total: every expected unit (an LD-score
## chromosome, a stratified GRM set, one REML fit) must be completed, QC-failed
## with a recorded reason, or failed. A missing status file leaves its units
## unaccounted and reconcile() stops (AGENTS.md 9).
##
## Metrics are continuous recovery against the REALIZED simulated h2 (bias,
## RMSE, 95% coverage, boundary or out-of-interval rate, Spearman with a
## cluster-bootstrap CI) plus the REML estimation-failure rate. No class and
## no threshold (AGENTS.md 2.3, 3).
##
## Usage: Rscript _h/06_summarize.R --run-id <id>

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02b_greml_simulation_benchmark", "_h")),
                 "greml_functions.R"))

opts <- parse_v2_args(require = "run_id")
run_dir <- run_dir_for(opts$run_id)
run <- list(run_id = opts$run_id, dir = run_dir)
cfg <- load_greml_config(run_dir)
man <- read_run_manifest(run_dir)
smoke <- identical(man[["smoke_run"]], "TRUE")
arm <- man[["arm"]]
modes <- fread(file.path(run_dir, "config", "reml-modes.tsv"))
n_boot <- as.integer(cfg$summaries$bootstrap_n)

read_status <- function(pattern) {
    f <- list.files(file.path(run_dir, "status"), pattern = pattern, full.names = TRUE)
    if (!length(f)) return(data.table(unit = character(), status = character(),
                                      reason = character()))
    rbindlist(lapply(f, fread, colClasses = "character"), fill = TRUE)
}
read_results <- function(pattern) {
    f <- list.files(file.path(run_dir, "results", "reml"), pattern = pattern,
                    full.names = TRUE)
    rbindlist(lapply(f, fread), fill = TRUE)
}

flags <- function(dt) {
    dt[, `:=`(ld_model = arm,
              simulated_phenotypes_only = TRUE,
              absolute_pve_interpretation_allowed_for_observed_loci = FALSE,
              transfers_to_cohort = identical(arm, cfg$arm2$id))]
    dt
}

if (identical(arm, cfg$arm1$id)) {
    ## ================================================================ arm 1
    sizes <- as.integer(strsplit(man[["sample_sizes"]], ",")[[1]])
    k <- as.integer(man[["n_phenotypes"]])
    ld_t <- fread(file.path(run_dir, "config", "ldscore-tasks.tsv"))
    expected <- CJ(N = sizes, pheno = seq_len(k), reml_mode = modes$reml_mode)
    expected[, unit := sprintf("%d:pheno_%d:%s", N, pheno, reml_mode)]
    exp_units <- c(paste0("ldscore|", sprintf("%d:%d", ld_t$N, ld_t$chr)),
                   paste0("stratify|", sizes),
                   paste0("reml|", expected$unit))
    st <- rbind(
        read_status("^ldscore_task_")[, unit := paste0("ldscore|", unit)],
        read_status("^stratify_task_")[, unit := paste0("stratify|", unit)],
        read_status("^reml_task_")[, unit := paste0("reml|", unit)])
    reconcile(exp_units, completed = st[status == "completed", unit],
              qc_failed = st[status == "qc_failed", unit],
              failed = st[status == "failed", unit], run = run,
              allow_failures = smoke)

    est <- read_results("^task_[0-9]+\\.tsv$")
    truth <- rbindlist(lapply(sizes, function(N) {
        m <- fread(file.path(repo_root(), man[["inputs_root"]], sprintf("sim_%d_indiv", N),
                             "snp_phenotype_mapping.tsv"))
        data.table(N = N, phenotype_id = m$phenotype_id,
                   truth_h2 = m[[cfg$arm1$truth_column]],
                   target_h2 = m[[cfg$arm1$target_column]],
                   num_causal_snps = m$num_causal_snps)
    }))
    est <- merge(est, truth, by = c("N", "phenotype_id"), all.x = TRUE)
    if (anyNA(est$truth_h2)) stop("A fitted phenotype has no truth row")
    expected[, status := st[match(paste0("reml|", expected$unit), unit), status]]
    fail <- expected[, .(n_expected = .N,
                         n_estimation_failed = sum(status == "qc_failed")),
                     by = .(N, reml_mode)]
    met <- est[, recovery_metrics(h2_hat, h2_se, truth_h2, constrained[1]),
               by = .(N, reml_mode)]
    met <- merge(fail, met, by = c("N", "reml_mode"), all.x = TRUE)
    met[, estimation_failure_rate := n_estimation_failed / n_expected]
    sp <- est[, spearman_cluster_ci(h2_hat, truth_h2, phenotype_id, n_boot,
                                    seed_for(cfg$seeds$namespace, "ar1", N, reml_mode)),
              by = .(N, reml_mode)]
    setnames(met, "N", "cell"); setnames(sp, "N", "cell")
    met[, `:=`(cell_type = "simulated_n", architecture = "v1_generator", h2_nominal = NA_real_)]
    sp[, `:=`(cell_type = "simulated_n", architecture = "v1_generator")]
} else {
    ## ================================================================ arm 2
    loci <- fread(file.path(run_dir, "config", "locus-tasks.tsv"))
    grid <- fread(file.path(run_dir, "config", "scenario-grid.tsv"))
    expected <- CJ(locus_index = loci$locus_index, g = seq_len(nrow(grid)),
                   reml_mode = modes$reml_mode)
    expected[, `:=`(vmr_id = loci$vmr_id[match(locus_index, loci$locus_index)],
                    h2_nominal = grid$h2[g], architecture = grid$architecture[g],
                    replicate = grid$replicate[g])]
    expected[, unit := sprintf("%s|h2=%s|%s|rep%d|%s", vmr_id, format(h2_nominal),
                               architecture, replicate, reml_mode)]
    st <- read_status("^reml_task_")
    reconcile(expected$unit, completed = st[status == "completed", unit],
              qc_failed = st[status == "qc_failed", unit],
              failed = st[status == "failed", unit], run = run,
              allow_failures = smoke)
    expected[, `:=`(status = st$status[match(unit, st$unit)],
                    reason = st$reason[match(unit, st$unit)])]
    est <- read_results("^task_[0-9]+\\.tsv$")
    geo <- read_results("^loci_task_[0-9]+\\.tsv$")
    write_atomic(geo, file.path(run_dir, "summary", "locus-geometry.tsv"))
    est[, truth_h2 := realized_h2]
    ## Locus-level QC failures are not estimator failures; keep them out of the
    ## estimation-failure denominator and report them on their own.
    expected[, locus_qc_failed := !is.na(reason) & startsWith(reason, "locus_")]
    fail <- expected[, .(n_expected = .N,
                         n_locus_qc_failed = sum(locus_qc_failed, na.rm = TRUE),
                         n_estimation_failed = sum(status == "qc_failed" & !locus_qc_failed,
                                                   na.rm = TRUE)),
                     by = .(reml_mode, architecture, h2_nominal)]
    met <- est[, recovery_metrics(h2_hat, h2_se, truth_h2, constrained[1]),
               by = .(reml_mode, architecture, h2_nominal)]
    met <- merge(fail, met, by = c("reml_mode", "architecture", "h2_nominal"), all.x = TRUE)
    met[, estimation_failure_rate := n_estimation_failed /
            pmax(1, n_expected - n_locus_qc_failed)]
    met[, `:=`(cell = man[["region"]], cell_type = "region_design_n")]
    sp <- rbind(
        est[, spearman_cluster_ci(h2_hat, truth_h2, vmr_id, n_boot,
                                  seed_for(cfg$seeds$namespace, man[["region"]],
                                           architecture, reml_mode)),
            by = .(reml_mode, architecture)],
        est[, spearman_cluster_ci(h2_hat, truth_h2, vmr_id, n_boot,
                                  seed_for(cfg$seeds$namespace, man[["region"]],
                                           "all", reml_mode)),
            by = .(reml_mode)][, architecture := "all"],
        use.names = TRUE)
    sp[, `:=`(cell = man[["region"]], cell_type = "region_design_n")]
}

met <- flags(met); sp <- flags(sp); est <- flags(est)
met[, reml_role := modes$role[match(reml_mode, modes$reml_mode)]]
sp[, reml_role := modes$role[match(reml_mode, modes$reml_mode)]]
write_atomic(est, file.path(run_dir, "summary", "estimates.tsv"))
write_atomic(met, file.path(run_dir, "summary", "recovery-metrics.tsv"))
write_atomic(sp, file.path(run_dir, "summary", "spearman.tsv"))
write_atomic(expected[, .(unit, status)], file.path(run_dir, "summary", "unit-status.tsv"))
append_manifest(run, list(n_completed_fits = nrow(est),
                          summarized_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
message("[02b] summarized ", nrow(est), " fits")
