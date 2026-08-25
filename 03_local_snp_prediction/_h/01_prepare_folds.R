#!/usr/bin/env Rscript
#### 03_local_snp_prediction -- donor fold assignment ####
##
## Usage:
##   Rscript _h/01_prepare_folds.R --run-id lsp-AA-caudate-20260817 [--allow-unlocked]
##
## Assigns donors to outer folds, once, for the whole run. This is a separate
## step from fitting for one reason: every VMR must use the SAME donor
## partition. If each array task drew its own folds, the pooled out-of-fold
## metrics across VMRs would be computed on incompatible splits, and the
## per-donor prediction counts required by config/prediction.yml would be
## meaningless.
##
## Folds are stratified only by things known without the phenotype. Stratifying
## on the outcome is a form of leakage.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "03_local_snp_prediction"

opts <- parse_v2_args(require = "run_id")
allow_unlocked <- isTRUE(opts$allow_unlocked)

run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)
manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]
    if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mval("cohort"); region <- mval("region")

prediction <- load_config("prediction")
assert_locked(list(prediction = prediction), allow_unlocked = allow_unlocked)

#' Read a PI-lockable key, refusing a null unless this is a smoke run.
locked_value <- function(cfg, key, smoke_default) {
    v <- tryCatch(config_get(cfg, key), error = function(e) NULL)
    if (!is.null(v)) return(v)
    if (!allow_unlocked) {
        stop("config/prediction.yml key '", key, "' is null and this is a ",
             "production run. The PI must lock it before 03 runs ",
             "(AGENTS.md 12/14).", call. = FALSE)
    }
    warning("Using smoke-test default for unlocked key '", key, "': ",
            smoke_default, call. = FALSE)
    smoke_default
}

n_outer   <- as.integer(locked_value(prediction, "folds.outer", 5L))
n_repeats <- as.integer(locked_value(prediction, "folds.repeats", 1L))

## ------------------------------------------------------------------ donors
## The donor list comes from the upstream VMR catalog run, not from a fresh
## phenotype read: 03 must evaluate exactly the donors 02 estimated on, or the
## two endpoints are not comparable.
vmr_run <- mval("upstream_vmr_catalog_run_id")
if (is.na(vmr_run) || !nzchar(vmr_run)) {
    stop("Run manifest carries no upstream_vmr_catalog_run_id; rerun 00_new_run.R")
}
## Module 01 writes the run-wide donor list as `vmr/donors_plink.txt`, two
## columns FID/IID. There is no single `results/vmr_meth.phen` -- phenotypes are
## one file per VMR under `vmr/phenotypes/`. Donors are keyed FID::IID, the same
## key `load_observed_locus()` aligns genotype rows on, so a fold assignment
## made here joins to a locus matrix there without any positional assumption.
donor_f <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs", vmr_run,
                     "vmr", "donors_plink.txt")
if (!file.exists(donor_f)) {
    stop("Donor list not found for upstream run ", vmr_run, ": ", donor_f,
         "\n  Check the accepted 01 run ID recorded in the manifest.")
}
donor_dt <- fread(donor_f, header = FALSE, colClasses = "character")
if (ncol(donor_dt) < 2L) stop("donors_plink.txt must have FID and IID columns")
donors <- paste(donor_dt[[1]], donor_dt[[2]], sep = "::")
assert_no_dups(donors, "donors in the Module 01 donor list")
assert_expected_n(length(donors), cohort, region)

## --------------------------------------------------------------- partition
## seed_for() makes the partition a deterministic function of (run_id, region,
## repeat) -- rerunning this script reproduces the identical assignment without
## anyone having to remember to save an RNG state.
folds <- rbindlist(lapply(seq_len(n_repeats), function(r) {
    set.seed(seed_for(opts$run_id, region = region, repeat_i = r))
    ## Random permutation into n_outer near-equal blocks. No outcome is consulted.
    assignment <- rep(seq_len(n_outer), length.out = length(donors))
    data.table(repeat_i = r,
               donor = donors,
               outer_fold = sample(assignment))
}))

## Guard the property the whole design rests on: every donor is held out exactly
## once per repeat.
chk <- folds[, .N, by = .(repeat_i, donor)]
if (any(chk$N != 1L)) {
    stop("Fold assignment is not a partition: ", sum(chk$N != 1L),
         " donor-repeat pairs appear more than once")
}

write_atomic(folds, file.path(run_dir, "donor-folds.tsv"))
append_manifest(list(dir = run_dir), list(
    n_donors = length(donors),
    donor_checksum = donor_checksum(donors),
    n_outer_folds = n_outer,
    n_repeats = n_repeats,
    folds_prepared_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
))

message("[03] ", length(donors), " donors x ", n_repeats, " repeat(s) x ",
        n_outer, " outer folds")
