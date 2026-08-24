#!/usr/bin/env Rscript
#### 04_repeat_repressive_architecture -- open a run ####
##
## Usage:
##   Rscript _h/00_new_run.R --cohort AA --region caudate [--allow-unlocked]
##
## 04 is the module the manuscript's central claim rests on, so the gate here is
## the strictest in the tree: it requires an accepted 02 run for this cell, and
## (when the secondary predictor is requested) an accepted 03 run as well.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "04_repeat_repressive_architecture"
MODULE_TAG <- "rra"

opts <- parse_v2_args(require = c("cohort", "region"))
allow_unlocked <- isTRUE(opts$allow_unlocked)
want_prediction <- !identical(opts$with_prediction, "FALSE")

annot <- load_config("repeat_annotations")
thresholds <- load_config("thresholds")
assert_locked(list(repeat_annotations = annot, thresholds = thresholds),
              allow_unlocked = allow_unlocked)

## The primary predictor is the continuous standardized within-cell local SNP
## contribution score among eligible loci. It preserves the validated ordering
## while making no claim that score distances are exact PVE differences.
upstream_02 <- require_accepted_upstream(
    "02_local_genetic_variance", opts$cohort, opts$region,
    allow_unaccepted = allow_unlocked)

upstream_03 <- if (want_prediction) {
    require_accepted_upstream("03_local_snp_prediction", opts$cohort, opts$region,
                              allow_unaccepted = allow_unlocked)
} else NULL

run <- new_run(
    module = MODULE_TAG, cohort = opts$cohort, region = opts$region,
    module_root = file.path(repo_root(), MODULE),
    vmr_set_id = upstream_02$vmr_set_id %||% NA_character_,
    upstream = c(
        list(local_genetic_variance_run_id = upstream_02$run_id %||% NA_character_),
        if (!is.null(upstream_03))
            list(local_snp_prediction_run_id = upstream_03$run_id %||% NA_character_)
        else list()
    ),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        config_repeat_annotations_sha256 = attr(annot, "config_sha256"),
        genome_build = annot$genome_build,
        repeatmasker_version = annot$repeatmasker$version %||% NA_character_,
        secondary_predictor = if (want_prediction) "r2_pred_oof" else "none"
    )
)

dir.create(file.path(run$dir, "results"), showWarnings = FALSE)
message("[04] run ", run$run_id, " opened")
cat(run$run_id, "\n")
