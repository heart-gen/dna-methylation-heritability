#!/usr/bin/env Rscript
#### 05_cpg_meqtl_burden -- acceptance gate ####
##
## Usage:
##   Rscript _h/04_check_burden.R --run-id cmb-AA-caudate-20260823
##
## AGENTS.md 7.5 names three things that must hold before this module's result
## is usable: the denominators must be audited, genomic inflation must be
## resolved before the figure freeze, and tested CpGs must be reported apart
## from prepared-but-untested ones. Each is a criterion here.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

MODULE <- "05_cpg_meqtl_burden"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]; if (length(v) == 0) NA_character_ else v[1]
}
smoke <- identical(mval("smoke_run"), "TRUE")
res_dir <- file.path(run_dir, "results")
meqtl <- load_config("meqtl_parameters")

audit <- fread(file.path(res_dir, "meqtl-denominator-audit.tsv"))
burden <- fread(file.path(res_dir, "vmr-meqtl-burden.tsv"))
recon <- fread(file.path(run_dir, "task_reconciliation.tsv"))
rn <- function(k) as.integer(recon$n[recon$category == k])

infl_f <- file.path(res_dir, "qc", "genomic-inflation.tsv")
lambda <- if (file.exists(infl_f)) fread(infl_f)$lambda[1] else NA_real_
lambda_min <- config_get(meqtl, "genomic_inflation.min")
lambda_max <- config_get(meqtl, "genomic_inflation.max")

## The burden fraction must never exceed 1, and its denominator must be the
## TESTED CpG count. A fraction above 1 means a VMR counted more significant
## CpGs than it had tested -- the exact denominator defect REVISION_GUIDE R8
## records in the legacy module.
frac_col <- if ("proportion_cpgs_with_sig_meqtl" %in% names(burden)) {
    "proportion_cpgs_with_sig_meqtl"
} else NA_character_

criteria <- data.table(
    criterion = c("chromosome_reconciliation_complete",
                  "cpg_accounting_balances",
                  "tested_reported_separately",
                  "burden_fraction_in_unit_interval",
                  "genomic_inflation_reported",
                  "genomic_inflation_resolved",
                  "continuous_predictor_is_primary"),
    passed = c(
        rn("unaccounted") == 0L && rn("unexpected") == 0L,
        audit$n_unaccounted == 0L,
        file.exists(file.path(res_dir, "untested-cpgs.tsv")) ||
            audit$n_prepared_but_untested == 0L,
        is.na(frac_col) ||
            all(burden[[frac_col]] >= 0 & burden[[frac_col]] <= 1, na.rm = TRUE),
        is.finite(lambda),
        is.finite(lambda) && lambda >= lambda_min && lambda <= lambda_max,
        identical(meqtl$predictability_score_column,
                  "local_snp_contribution_score_z")
    ),
    detail = c(
        sprintf("unaccounted=%d unexpected=%d", rn("unaccounted"), rn("unexpected")),
        sprintf("prepared=%d tested=%d untested=%d unaccounted=%d",
                audit$n_prepared_cpgs, audit$n_tested_cpgs,
                audit$n_prepared_but_untested, audit$n_unaccounted),
        sprintf("prepared_but_untested=%d", audit$n_prepared_but_untested),
        if (is.na(frac_col)) "no burden fraction column" else
            sprintf("range=[%.3f, %.3f]", min(burden[[frac_col]], na.rm = TRUE),
                    max(burden[[frac_col]], na.rm = TRUE)),
        sprintf("lambda=%s", format(lambda, digits = 4)),
        sprintf("lambda=%s (allowed %s-%s)", format(lambda, digits = 4),
                lambda_min, lambda_max),
        meqtl$predictability_score_column
    )
)

decision <- if (all(criteria$passed)) {
    if (smoke) "PASS_SMOKE_ONLY_NOT_ACCEPTABLE" else "PASS_CPG_MEQTL_BURDEN_QC"
} else "FAIL_CPG_MEQTL_BURDEN_QC"

criteria[, `:=`(run_id = opts$run_id, region = mval("region"),
                population = mval("cohort"), decision = decision)]
write_atomic(criteria, file.path(res_dir, "burden-qc-criteria.tsv"))
write_atomic(data.table(run_id = opts$run_id, decision = decision,
                        smoke_run = smoke, n_criteria = nrow(criteria),
                        n_passed = sum(criteria$passed)),
             file.path(res_dir, "burden-decision.tsv"))
print(criteria[, .(criterion, passed, detail)])
message("[05] decision: ", decision)
if (!all(criteria$passed)) quit(status = 1L)
