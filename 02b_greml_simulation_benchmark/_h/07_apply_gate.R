#!/usr/bin/env Rscript
#### 02b -- acceptance gate ####
##
## PASS_GREML_BENCHMARK_QC certifies completeness and provenance, NOT good
## recovery. Poor recovery at small n is a legitimate result -- it is the result
## this module exists to measure -- and must not block sealing, as Module 06's
## null and Module 10's coverage gate did not.
##
## Usage: Rscript _h/07_apply_gate.R --run-id <id>

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
arm1 <- identical(man[["arm"]], cfg$arm1$id)

rec <- fread(file.path(run_dir, "task_reconciliation.tsv"))
n_of <- function(cat) rec[category == cat, n]
est <- fread(file.path(run_dir, "summary", "estimates.tsv"))
met <- fread(file.path(run_dir, "summary", "recovery-metrics.tsv"))
primary <- fread(file.path(run_dir, "config", "reml-modes.tsv"))[role == "primary", reml_mode]
ep <- est[reml_mode == primary]

inputs_pinned <- if (arm1) {
    f <- file.path(run_dir, "inputs", "input_checksums.tsv")
    file.exists(f) && nrow(fread(f)) ==
        5L * length(strsplit(man[["sample_sizes"]], ",")[[1]])
} else {
    up <- tryCatch(require_accepted_upstream(cfg$arm2$upstream, man[["cohort"]],
                                             man[["region"]]),
                   error = function(e) list(run_id = NA_character_))
    identical(up$run_id, man[["upstream_vmr_catalog_run_id"]])
}
truth_covered <- if (arm1) {
    nrow(ep[truth_h2 < 0.01]) > 0 && nrow(ep[truth_h2 > 0.5]) > 0
} else {
    h2v <- as.numeric(unlist(cfg$arm2$h2_values))
    nrow(ep[h2_nominal == min(h2v)]) > 0 && nrow(ep[h2_nominal == max(h2v)]) > 0
}
crit <- data.table(
    criterion = c("units_reconciled_zero_failures",
                  "gcta_version_matches_pin",
                  if (arm1) "simulated_inputs_checksummed" else "upstream_catalog_still_accepted",
                  "truth_range_covered_by_primary_fits",
                  "simulated_only_flags_on_every_row",
                  "primary_mode_metrics_and_failure_rate_reported"),
    passed = c(
        n_of("failed") == 0 && n_of("unaccounted") == 0 && n_of("unexpected") == 0,
        identical(man[["gcta_version"]], as.character(cfg$gcta$version)) &&
            identical(gcta_version(cfg$gcta$binary), as.character(cfg$gcta$version)),
        isTRUE(inputs_pinned),
        isTRUE(truth_covered),
        nrow(est) > 0 && all(est$simulated_phenotypes_only %in% TRUE) &&
            all(est$absolute_pve_interpretation_allowed_for_observed_loci %in% FALSE),
        nrow(met[reml_mode == primary]) > 0 &&
            !anyNA(met[reml_mode == primary, estimation_failure_rate])))
decision <- if (all(crit$passed)) "PASS_GREML_BENCHMARK_QC" else "FAIL_GREML_BENCHMARK_QC"
crit[, `:=`(run_id = opts$run_id, decision = decision, smoke_run = smoke)]
print(crit)
write_atomic(crit, file.path(run_dir, "summary", "gate-criteria.tsv"))
append_manifest(run, list(greml_benchmark_decision = decision,
                          gate_criteria_passed = sprintf("%d/%d", sum(crit$passed), nrow(crit))))
