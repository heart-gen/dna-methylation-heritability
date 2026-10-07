#!/usr/bin/env Rscript
#### 05_cpg_meqtl_burden -- open the cross-region slope-inference run ####
##
## Usage:
##   Rscript _h/06_slope_inference_new_run.R --cohort AA [--allow-unlocked]
##
## Opens `cmb-{cohort}-crossregion-{date}` and writes the paired donor-bootstrap
## draws every later stage reads. Why the run exists is in
## 08_slope_inference.R; this stage only fixes WHICH donors each draw holds, so
## that the draws are decided once, before any mapping, and can be audited.
##
## The bootstrap is PAIRED: one draw resamples the UNION of the contrast
## regions' donors, and each region then maps the drawn donors it has, with
## their multiplicity. DLPFC and hippocampus share 115 of 118 donors, so their
## slopes are strongly positively correlated, and only a paired draw measures
## that correlation instead of assuming it away -- the construction
## 09b_aging_application/_h/05_cross_region_concordance.R uses for the same
## contrast. Resampling is stratified by diagnosis, as there, because cases and
## controls are not exchangeable in this cohort.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages(library(data.table))

MODULE <- "05_cpg_meqtl_burden"
MODULE_TAG <- "cmb"

opts <- parse_v2_args(require = "cohort")
allow_unlocked <- isTRUE(opts$allow_unlocked)
cohort <- opts$cohort

meqtl <- load_config("meqtl_parameters")
si <- meqtl$cross_region_slope_inference
if (is.null(si)) stop("config/meqtl_parameters.yml has no cross_region_slope_inference block")
contrast <- as.character(unlist(si$contrast))
if (length(contrast) != 2L) stop("cross_region_slope_inference.contrast must name two regions")
confounded <- as.character(config_get(load_config("region_donor_generalization"),
    "interpretation.technically_confounded_regions.all_outcomes"))
if (length(intersect(contrast, confounded))) {
    stop("The contrast includes a technically confounded region; tier 2 exists ",
         "because it avoids the batch boundary (AGENTS.md 7.7, 8.1).")
}
B <- as.integer(if (allow_unlocked) si$smoke_bootstrap_n else si$bootstrap_n)

## ------------------------------------------------------------ the two cells
up <- lapply(contrast, function(re)
    require_accepted_upstream(MODULE, cohort, re, allow_unaccepted = allow_unlocked))
names(up) <- contrast
if (anyNA(vapply(up, function(u) as.character(u$run_id), character(1)))) {
    stop("Slope inference reads the accepted cells' prepared inputs; there is ",
         "no smoke input path.")
}
cell_dir <- function(re) file.path(repo_root(), MODULE, "_m", "runs", up[[re]]$run_id)

## The stratifier is a covariate the cells actually fitted, read from their own
## prepared design (rows = covariates, columns = donors). Every autosome carries
## the same donor set; chr1 is read and the others are checked against it in the
## mapping stage, which loads them anyway.
strat_col <- si$stratify_by
donors <- rbindlist(lapply(contrast, function(re) {
    cov <- fread(file.path(cell_dir(re), "inputs", "chr1.covariates.tsv"))
    setnames(cov, 1, "covariate")
    row <- cov[covariate == strat_col]
    if (nrow(row) != 1L) stop(up[[re]]$run_id, " has no covariate row '", strat_col, "'")
    ids <- setdiff(names(cov), "covariate")
    data.table(region = re, donor = ids, stratum = as.character(unlist(row[, ..ids])))
}))
union <- unique(donors[, .(donor, stratum)])
if (anyDuplicated(union$donor)) {
    stop("A donor carries different '", strat_col, "' values in the two regions")
}
setorder(union, donor)

## -------------------------------------------------------------- the draws
seed <- as.integer(si$seed)
draws <- rbindlist(lapply(seq_len(B), function(b) {
    set.seed(seed + b)
    picked <- unlist(lapply(split(union$donor, union$stratum), function(ids)
        ids[sample.int(length(ids), length(ids), replace = TRUE)]), use.names = FALSE)
    data.table(draw = b, slot = seq_along(picked), donor = picked)
}))

if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked)")
}
run <- new_run(
    module = MODULE_TAG, cohort = cohort, region = "crossregion",
    module_root = file.path(repo_root(), MODULE), run_id = opts$run_id,
    upstream = stats::setNames(lapply(contrast, function(re) up[[re]]$run_id),
                               paste0("cpg_meqtl_burden_", contrast)),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        config_meqtl_sha256 = attr(meqtl, "config_sha256"),
        analysis_kind = "cross_region_slope_inference",
        contrast = paste(contrast, collapse = ","),
        bootstrap_n = B,
        bootstrap_nperm = si$bootstrap_nperm,
        bootstrap_seed = seed,
        stratify_by = strat_col,
        n_union_donors = nrow(union),
        n_shared_donors = length(Reduce(intersect, split(donors$donor, donors$region)))
    )
)
dir.create(file.path(run$dir, "inputs"), showWarnings = FALSE)
dir.create(file.path(run$dir, "results", "bootstrap"), recursive = TRUE, showWarnings = FALSE)
write_atomic(union, file.path(run$dir, "inputs", "bootstrap-donors.tsv"))
write_atomic(draws, file.path(run$dir, "inputs", "bootstrap-draws.tsv"))
message("[05] slope-inference run ", run$run_id, ": ", B, " paired draws over ",
        nrow(union), " donors")
cat(run$run_id, "\n", sep = "")
