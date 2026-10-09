#!/usr/bin/env Rscript
#### 05_cpg_meqtl_burden -- open the cross-region slope-inference run ####
##
## Usage:
##   Rscript _h/06_slope_inference_new_run.R --cohort AA [--allow-unlocked]
##
## Opens `cmb-{cohort}-crossregion-{date}` and writes the paired delete-d
## draws every later stage reads. Why the run exists is in
## 08_slope_inference.R; this stage only fixes WHICH donors each draw deletes,
## so that the draws are decided once, before any mapping, and can be audited.
##
## The draws are PAIRED: one draw deletes d donors from the UNION of the
## contrast regions' donors, and each region then maps the donors it has left.
## DLPFC and hippocampus share 115 of 118 donors, so their slopes are strongly
## positively correlated, and only a paired draw measures that correlation
## instead of assuming it away.
##
## Donors are drawn WITHOUT replacement. A bootstrap was built first and its
## smoke run called every tested CpG significant in every draw: a duplicated
## donor is duplicated in genotype and phenotype alike, the permutation null
## shuffles the duplicates apart, and the null comes out too narrow
## (config/meqtl_parameters.yml records the numbers). A subset of distinct
## donors keeps every permutation valid.

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
B <- as.integer(if (allow_unlocked) si$smoke_draws_n else si$draws_n)
if (!identical(si$resampling, "paired_delete_d_donor_jackknife")) {
    stop("cross_region_slope_inference.resampling must be paired_delete_d_donor_jackknife")
}

## ------------------------------------------------------------ the two cells
up <- lapply(contrast, function(re)
    require_accepted_upstream(MODULE, cohort, re, allow_unaccepted = allow_unlocked))
names(up) <- contrast
if (anyNA(vapply(up, function(u) as.character(u$run_id), character(1)))) {
    stop("Slope inference reads the accepted cells' prepared inputs; there is ",
         "no smoke input path.")
}
cell_dir <- function(re) file.path(repo_root(), MODULE, "_m", "runs", up[[re]]$run_id)

## Donor sets come from each cell's own prepared design (rows = covariates,
## columns = donors). Every autosome carries the same donor set; chr1 is read
## here and 07_subsample_map.py intersects each chromosome's genotypes anyway.
donors <- rbindlist(lapply(contrast, function(re) {
    cov <- fread(file.path(cell_dir(re), "inputs", "chr1.covariates.tsv"), nrows = 1L)
    data.table(region = re, donor = setdiff(names(cov), names(cov)[1]))
}))
union <- data.table(donor = sort(unique(donors$donor)))
union[, (paste0("in_", contrast)) := lapply(contrast, function(re)
    donor %in% donors[region == re, donor])]

## -------------------------------------------------------------- the draws
n_u <- nrow(union)
d <- as.integer(round(as.numeric(si$delete_fraction) * n_u))
if (d <= sqrt(n_u) || d >= n_u / 2) {
    stop("delete_fraction gives d = ", d, " of n = ", n_u, "; the delete-d ",
         "jackknife needs sqrt(n) < d, and d < n/2 keeps each draw's mapping ",
         "power near the full sample's")
}
seed <- as.integer(si$seed)
draws <- rbindlist(lapply(seq_len(B), function(b) {
    set.seed(seed + b)
    data.table(draw = b, donor = sort(union$donor[sample.int(n_u, d)]))
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
        resampling = si$resampling,
        draws_n = B,
        delete_d = d,
        draw_nperm = si$draw_nperm,
        draw_seed = seed,
        n_union_donors = n_u,
        n_shared_donors = length(Reduce(intersect, split(donors$donor, donors$region)))
    )
)
dir.create(file.path(run$dir, "inputs"), showWarnings = FALSE)
dir.create(file.path(run$dir, "results", "subsample"), recursive = TRUE, showWarnings = FALSE)
write_atomic(union, file.path(run$dir, "inputs", "subsample-donors.tsv"))
## One row per DELETED donor: a draw's mapped set is the union minus these.
write_atomic(draws, file.path(run$dir, "inputs", "subsample-deletions.tsv"))
message("[05] slope-inference run ", run$run_id, ": ", B, " paired draws, each ",
        "deleting ", d, " of ", n_u, " donors")
cat(run$run_id, "\n", sep = "")
