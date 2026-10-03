#!/usr/bin/env Rscript
#### 02c -- acceptance gate ####
##
## PASS_CIS_GREML_SENSITIVITY_QC certifies completeness, provenance and that
## 02c fitted the same loci and SNPs as Module 02. It is not a success
## criterion: a weak or null ordering is a legitimate result and does not
## block sealing, and nothing about the Module 02 score changes either way.
##
## Usage: Rscript _h/03_apply_gate.R --run-id <id>

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02c_cis_greml_sensitivity", "_h")),
                 "cgs_functions.R"))

opts <- parse_v2_args(require = "run_id")
run_dir <- run_dir_for(opts$run_id)
run <- list(run_id = opts$run_id, dir = run_dir)
cfg <- load_cgs_config(run_dir)
man <- read_run_manifest(run_dir)
smoke <- identical(man[["smoke_run"]], "TRUE")
rec <- fread(file.path(run_dir, "task_reconciliation.tsv"))
n_of <- function(x) rec[category == x, n]
loci <- fread(file.path(run_dir, "summary", "locus-snp-check.tsv"))
vmr <- fread(file.path(run_dir, "summary", "vmr-cis-greml.tsv"))
ex <- fread(file.path(run_dir, "summary", "existence.tsv"))
od <- fread(file.path(run_dir, "summary", "ordering.tsv"))
up_ok <- all(vapply(names(cfg$upstreams), function(k) {
    u <- tryCatch(require_accepted_upstream(cfg$upstreams[[k]], man[["cohort"]],
                                            man[["region"]]), error = function(e) NULL)
    !is.null(u) && identical(u$run_id, man[[paste0("upstream_", k, "_run_id")]])
}, logical(1)))
chk <- loci[locus_status == "ok" & !is.na(snp_count_matches_module02)]
crit <- data.table(
    criterion = c("units_reconciled_zero_failures", "gcta_version_matches_pin",
                  "upstreams_01_and_02_still_accepted",
                  "snp_set_matches_module02_for_every_fitted_vmr",
                  "interpretation_flags_on_every_row",
                  "primary_existence_and_ordering_reported"),
    passed = c(
        n_of("failed") == 0 && n_of("unaccounted") == 0 && n_of("unexpected") == 0,
        identical(man[["gcta_version"]], as.character(cfg$gcta$version)) &&
            identical(gcta_version(cfg$gcta$binary), as.character(cfg$gcta$version)),
        up_ok,
        nrow(chk) > 0 && all(chk$snp_count_matches_module02 %in% TRUE),
        nrow(vmr) > 0 && all(vmr$replaces_module_02_score %in% FALSE) &&
            all(vmr$per_vmr_absolute_h2_reportable %in% FALSE) &&
            all(vmr$greml_significance_class_allowed %in% FALSE),
        is.finite(ex[reml_role == "primary", estimate]) &&
            is.finite(od[reml_role == "primary" & role_of_comparison == "primary_ordering",
                         estimate])))
decision <- if (all(crit$passed)) "PASS_CIS_GREML_SENSITIVITY_QC" else
    "FAIL_CIS_GREML_SENSITIVITY_QC"
crit[, `:=`(run_id = opts$run_id, decision = decision, smoke_run = smoke,
            n_snp_set_checked = nrow(chk),
            n_snp_set_mismatch = sum(!(chk$snp_count_matches_module02 %in% TRUE)))]
print(crit)
write_atomic(crit, file.path(run_dir, "summary", "gate-criteria.tsv"))
append_manifest(run, list(cis_greml_sensitivity_decision = decision,
                          gate_criteria_passed = sprintf("%d/%d", sum(crit$passed), nrow(crit))))
