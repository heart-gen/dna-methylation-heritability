#!/usr/bin/env Rscript
#### 07_transcription_splicing_coupling -- open a run ####
##
## Usage:
##   Rscript _h/00_new_run.R --cohort AA --region caudate [--allow-unlocked]
##
## Three upstreams are consumed directly: 01 for the corrected VMR boundaries
## and per-donor VMR methylation, 02 for the relative local-control score, and
## 05 for the CpG meQTL burden. All three must describe the same vmr_set_id, or
## the coupling test would join methylation, score and burden for different loci
## that happen to share an ID.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "07_transcription_splicing_coupling"
MODULE_TAG <- "tsc"

opts <- parse_v2_args(require = c("cohort", "region"))
allow_unlocked <- isTRUE(opts$allow_unlocked)

ts <- load_config("transcription_splicing")
assert_locked(list(transcription_splicing = ts), allow_unlocked = allow_unlocked)

if (!identical(ts$coupling$predictors$local_genetic_control,
               "local_snp_contribution_score_z")) {
    stop("config/transcription_splicing.yml must use ",
         "local_snp_contribution_score_z as the local-genetic-control ",
         "predictor; legacy predictability and absolute h2 are banned ",
         "(AGENTS.md 3).")
}

## AGENTS.md 7.6 forbids an unbounded transcriptome-wide screen. The tested
## universe is the enabled modality list, and it is frozen here so a modality
## cannot be added after results are seen.
enabled <- names(ts$modalities)[vapply(ts$modalities,
                                       function(m) isTRUE(m$enabled), logical(1))]
if (length(enabled) == 0) stop("No modality is enabled in config")

## Assay files must exist for this region before a run directory is created.
for (mod in enabled) {
    assay <- ts$modalities[[mod]]$assay
    f <- file.path(repo_root(), ts$assay_files[[assay]][[opts$region]])
    if (!file.exists(f)) {
        stop("Missing ", assay, " assay for ", opts$region, ": ", f)
    }
}

upstream_01 <- require_accepted_upstream(
    "01_vmr_catalog", opts$cohort, opts$region,
    allow_unaccepted = allow_unlocked)
upstream_02 <- require_accepted_upstream(
    "02_local_genetic_variance", opts$cohort, opts$region,
    allow_unaccepted = allow_unlocked)
upstream_05 <- require_accepted_upstream(
    "05_cpg_meqtl_burden", opts$cohort, opts$region,
    allow_unaccepted = allow_unlocked)

sets <- c(upstream_01$vmr_set_id, upstream_02$vmr_set_id, upstream_05$vmr_set_id)
sets <- sets[!is.na(sets)]
if (length(unique(sets)) > 1) {
    stop("vmr_set_id mismatch across upstreams: ",
         paste(unique(sets), collapse = " vs "),
         "\n  01, 02 and 05 must describe the same VMR set (AGENTS.md 6).")
}

if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked); ",
         "production run IDs are derived, not chosen.")
}

run <- new_run(
    module = MODULE_TAG, cohort = opts$cohort, region = opts$region,
    module_root = file.path(repo_root(), MODULE),
    run_id = opts$run_id,
    vmr_set_id = if (length(sets)) sets[1] else NA_character_,
    upstream = list(
        vmr_catalog_run_id = upstream_01$run_id %||% NA_character_,
        local_genetic_variance_run_id = upstream_02$run_id %||% NA_character_,
        cpg_meqtl_burden_run_id = upstream_05$run_id %||% NA_character_
    ),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        config_transcription_splicing_sha256 = attr(ts, "config_sha256"),
        modalities = paste(enabled, collapse = ","),
        fdr_family = ts$association$fdr_family,
        fdr_threshold = ts$association$fdr_threshold,
        internal_libd_eqtl_support_arm =
            if (isTRUE(ts$internal_libd_eqtl_support_arm)) "TRUE" else "FALSE"
    )
)

dir.create(file.path(run$dir, "results"), showWarnings = FALSE)
dir.create(file.path(run$dir, "links"), showWarnings = FALSE)
message("[07] run ", run$run_id, " opened (modalities: ",
        paste(enabled, collapse = ", "), ")")
cat(run$run_id, "\n", sep = "")
