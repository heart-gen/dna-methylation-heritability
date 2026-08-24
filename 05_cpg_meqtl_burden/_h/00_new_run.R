#!/usr/bin/env Rscript
#### 05_cpg_meqtl_burden -- open a run ####
##
## Usage:
##   Rscript _h/00_new_run.R --cohort AA --region caudate [--allow-unlocked]
##
## 05 consumes BOTH upstreams directly: 02 for the relative local-control
## gradient, and 01 for the corrected CpG-to-VMR membership table. The membership
## table is the part that matters most -- the legacy burden analysis aggregated
## CpGs to VMRs using a membership built by the pre-V1-fix pipeline, so its
## denominators were wrong independently of the meQTL mapping itself.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "05_cpg_meqtl_burden"
MODULE_TAG <- "cmb"

opts <- parse_v2_args(require = c("cohort", "region"))
allow_unlocked <- isTRUE(opts$allow_unlocked)

meqtl <- load_config("meqtl_parameters")
thresholds <- load_config("thresholds")
assert_locked(list(thresholds = thresholds), allow_unlocked = allow_unlocked)

## The primary predictor is the validated relative local-control score. Refuse
## older absolute-h2 or legacy predictability configurations.
if (!identical(meqtl$predictability_score_column,
               "local_snp_contribution_score_z")) {
    stop("config/meqtl_parameters.yml sets predictability_score_column to '",
         meqtl$predictability_score_column,
         "'. v2 requires local_snp_contribution_score_z; absolute-h2 and ",
         "legacy predictability inputs are prohibited.")
}

upstream_02 <- require_accepted_upstream(
    "02_local_genetic_variance", opts$cohort, opts$region,
    allow_unaccepted = allow_unlocked)
upstream_01 <- require_accepted_upstream(
    "01_vmr_catalog", opts$cohort, opts$region,
    allow_unaccepted = allow_unlocked)

## The two upstreams must describe the SAME VMR set, or CpG membership and score
## estimates refer to different loci with colliding IDs.
if (!is.na(upstream_02$vmr_set_id) && !is.na(upstream_01$vmr_set_id) &&
    !identical(upstream_02$vmr_set_id, upstream_01$vmr_set_id)) {
    stop("vmr_set_id mismatch between upstreams: 02 cites ",
         upstream_02$vmr_set_id, ", 01 cites ", upstream_01$vmr_set_id)
}

run <- new_run(
    module = MODULE_TAG, cohort = opts$cohort, region = opts$region,
    module_root = file.path(repo_root(), MODULE),
    vmr_set_id = upstream_01$vmr_set_id %||% NA_character_,
    upstream = list(
        vmr_catalog_run_id = upstream_01$run_id %||% NA_character_,
        local_genetic_variance_run_id = upstream_02$run_id %||% NA_character_
    ),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        config_meqtl_sha256 = attr(meqtl, "config_sha256"),
        cis_window_bp = meqtl$cis_window_bp,
        mapping_engine = meqtl$mapping$engine,
        fdr_threshold = meqtl$mapping$fdr_threshold
    )
)

dir.create(file.path(run$dir, "results"), showWarnings = FALSE)
message("[05] run ", run$run_id, " opened")
cat(run$run_id, "\n")
