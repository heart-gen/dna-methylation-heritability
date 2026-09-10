#!/usr/bin/env Rscript
#### 03_local_snp_prediction -- open a run ####
##
## Usage:
##   Rscript _h/00_new_run.R --cohort AA --region caudate [--allow-unlocked]
##
## Creates _m/runs/{RUN_ID}/, records the upstream 02 run it will consume, and
## writes the VMR task manifest the array steps index into. Nothing scientific
## happens here; this exists so that the gate is checked ONCE, up front, rather
## than 11,000 times inside array tasks that have already been queued.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "03_local_snp_prediction"
MODULE_TAG <- "lsp"

opts <- parse_v2_args(require = c("cohort", "region"))
allow_unlocked <- isTRUE(opts$allow_unlocked)

prediction <- load_config("prediction")
thresholds <- load_config("thresholds")
assert_locked(list(prediction = prediction, thresholds = thresholds),
              allow_unlocked = allow_unlocked)

## ---------------------------------------------------------------- gate 02
## AGENTS.md 6. The gate remains closed until an observed relative-score run is
## accepted. That is intended behavior, not a bug to work around.
upstream <- require_accepted_upstream(
    "02_local_genetic_variance", opts$cohort, opts$region,
    allow_unaccepted = allow_unlocked)

lcg <- load_local_genetic_control(
    upstream$run_id, region = opts$region, cohort = opts$cohort,
    eligible_only = FALSE
)

## Prediction is evaluated on every locus 02 produced a summary for, including
## loci ineligible for relative-score interpretation: predictive accuracy is an
## empirical held-out quantity and does not depend on the calibration domain.
tasks <- data.table(
    task_id = seq_len(nrow(lcg)),
    vmr_id  = lcg$vmr_id,
    chrom   = lcg$chrom,
    start   = lcg$start,
    end     = lcg$end,
    local_genetic_control_eligible = lcg$local_genetic_control_eligible,
    local_genetic_control_exclusion_reason =
        lcg$local_genetic_control_exclusion_reason
)
assert_no_dups(tasks$vmr_id, "VMR IDs from the upstream 02 table")

## A smoke run must be nameable, so it can never be mistaken for the production
## run of the same cell, and sizeable, so the chain can be exercised in minutes.
## Both are refused unless --allow-unlocked marks this a smoke run.
if (!is.null(opts$smoke_n)) {
    if (!allow_unlocked) stop("--smoke-n requires --allow-unlocked")
    keep <- unique(round(seq(1, nrow(tasks), length.out = as.integer(opts$smoke_n))))
    tasks <- tasks[keep]
    tasks[, task_id := seq_len(.N)]
}
if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked); ",
         "production run IDs are derived, not chosen.")
}

## ------------------------------------------------- cell identity from 02
## 03 must score the SAME donors and the SAME variants as 02 (that is why both
## call one load_observed_locus()). The cell identity therefore travels from
## 02's manifest rather than being re-derived here, where it could drift.
## Pre-2026-09-10 runs of 02 carry none of these fields; the fallbacks are the
## behaviour those runs actually had.
lgv_manifest_path <- file.path(repo_root(), "02_local_genetic_variance",
                               "_m", "runs", upstream$run_id, "manifest.tsv")
lgv_field <- function(field, default = NA_character_) {
    if (!file.exists(lgv_manifest_path)) return(default)
    m <- fread(lgv_manifest_path, colClasses = "character")
    v <- m$value[m$field == field]
    if (length(v) == 1L && !is.na(v) && nzchar(v)) as.character(v[[1L]]) else default
}
catalog_cohort   <- lgv_field("catalog_cohort", opts$cohort)
estimation_group <- lgv_field("estimation_group", opts$cohort)
covar_prefix     <- lgv_field("covar_prefix",
                              if (identical(opts$cohort, "AA"))
                                  "TOPMed_LIBD.AA" else "TOPMed_LIBD")
upstream_module  <- lgv_field("upstream_module", "01_vmr_catalog")

run <- new_run(
    module = MODULE_TAG, cohort = opts$cohort, region = opts$region,
    module_root = file.path(repo_root(), MODULE),
    run_id = opts$run_id,
    vmr_set_id = upstream$vmr_set_id %||% NA_character_,
    upstream = list(
        local_genetic_variance_run_id = upstream$run_id %||% NA_character_,
        vmr_catalog_run_id = attr(lcg, "upstream_vmr_run_id") %||%
            (if ("upstream_vmr_run_id" %in% names(lcg)) lcg$upstream_vmr_run_id[1] else NA_character_)
    ),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        catalog_cohort = catalog_cohort,
        estimation_group = estimation_group,
        covar_prefix = covar_prefix,
        upstream_vmr_module = upstream_module,
        config_prediction_sha256 = attr(prediction, "config_sha256"),
        evaluation_standard = prediction$evaluation_standard,
        n_expected_tasks = nrow(tasks)
    )
)

dir.create(file.path(run$dir, "results"), showWarnings = FALSE)
write_atomic(tasks, file.path(run$dir, "task-manifest.tsv"))

message("[03] run ", run$run_id, " opened with ", nrow(tasks), " VMR tasks")
cat(run$run_id, "\n", sep = "")   # sep="": cat() otherwise emits "<id> \n" and the shell captures the trailing space
