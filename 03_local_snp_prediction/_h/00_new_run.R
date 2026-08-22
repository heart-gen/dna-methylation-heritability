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

run <- new_run(
    module = MODULE_TAG, cohort = opts$cohort, region = opts$region,
    module_root = file.path(repo_root(), MODULE),
    vmr_set_id = upstream$vmr_set_id %||% NA_character_,
    upstream = list(
        local_genetic_variance_run_id = upstream$run_id %||% NA_character_,
        vmr_catalog_run_id = attr(lcg, "upstream_vmr_run_id") %||%
            (if ("upstream_vmr_run_id" %in% names(lcg)) lcg$upstream_vmr_run_id[1] else NA_character_)
    ),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        config_prediction_sha256 = attr(prediction, "config_sha256"),
        evaluation_standard = prediction$evaluation_standard,
        n_expected_tasks = nrow(tasks)
    )
)

dir.create(file.path(run$dir, "results"), showWarnings = FALSE)
write_atomic(tasks, file.path(run$dir, "task-manifest.tsv"))

message("[03] run ", run$run_id, " opened with ", nrow(tasks), " VMR tasks")
cat(run$run_id, "\n")
