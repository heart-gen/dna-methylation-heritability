#!/usr/bin/env Rscript
#### 06_partitioned_heritability -- open a run ####
##
## Usage:
##   Rscript _h/00_new_run.R --cohort AA --region caudate [--allow-unlocked]
##
## AGENTS.md 6 carves this module out of the 03-05 chain explicitly: "06_
## partitioned_heritability depends only on an accepted 02_local_genetic_variance
## score; its position in this list is a total order, not a claim that it
## consumes 03-05." So exactly one upstream gate is checked here, and checking
## more would be wrong rather than merely cautious.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "06_partitioned_heritability"
MODULE_TAG <- "sldsc"

opts <- parse_v2_args(require = c("cohort", "region"))
allow_unlocked <- isTRUE(opts$allow_unlocked)

ph <- load_config("partitioned_heritability")
assert_locked(list(partitioned_heritability = ph), allow_unlocked = allow_unlocked)

## The annotation may be built from one column only. A config that has drifted
## back to an absolute-PVE or legacy-predictability column is a scientific
## error, not a preference, so it stops the run before a run directory exists.
if (!identical(ph$annotation$score_column, "local_snp_contribution_score_z")) {
    stop("config/partitioned_heritability.yml sets annotation.score_column to '",
         ph$annotation$score_column,
         "'. Module 06 requires local_snp_contribution_score_z (AGENTS.md 3).")
}
if (!isTRUE(ph$annotation$continuous) ||
    !isTRUE(ph$annotation$forbid_thresholding) ||
    !isTRUE(ph$annotation$forbid_grouping)) {
    stop("The annotation must be continuous, unthresholded and ungrouped. ",
         "The retired v1 quintile form is banned by AGENTS.md 3.")
}

arm <- ph$ld_reference_arm
if (is.null(ph$ld_references[[arm]])) {
    stop("ld_reference_arm '", arm, "' has no entry under ld_references. ",
         "The AFR sensitivity panel is built by a separate issue and is not ",
         "populated yet.")
}

## The frozen trait list defines the FDR family (AGENTS.md 7: the tested
## universe is declared before results). Refuse to open a run whose declared
## traits are not all present on disk -- a silently skipped trait would shrink
## the family after the fact and inflate every remaining q-value.
missing_traits <- vapply(ph$traits, function(t) !file.exists(t$file), logical(1))
if (any(missing_traits)) {
    stop("Declared GWAS summary statistics missing for: ",
         paste(vapply(ph$traits[missing_traits], `[[`, character(1), "name"),
               collapse = ", "),
         "\n  The trait list is the FDR family and may not be quietly reduced.")
}

upstream_02 <- require_accepted_upstream(
    "02_local_genetic_variance", opts$cohort, opts$region,
    allow_unaccepted = allow_unlocked)

if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked); ",
         "production run IDs are derived, not chosen.")
}

run <- new_run(
    module = MODULE_TAG, cohort = opts$cohort, region = opts$region,
    module_root = file.path(repo_root(), MODULE),
    run_id = opts$run_id,
    vmr_set_id = upstream_02$vmr_set_id %||% NA_character_,
    upstream = list(
        local_genetic_variance_run_id = upstream_02$run_id %||% NA_character_
    ),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        config_partitioned_heritability_sha256 = attr(ph, "config_sha256"),
        annotation_score_column = ph$annotation$score_column,
        ld_reference_arm = arm,
        ld_reference_label = ph$ld_references[[arm]]$label,
        n_traits = length(ph$traits),
        fdr_family = ph$fdr_family,
        fdr_method = ph$fdr_method
    )
)

dir.create(file.path(run$dir, "results"), showWarnings = FALSE)
dir.create(file.path(run$dir, "annotation"), showWarnings = FALSE)
dir.create(file.path(run$dir, "sumstats"), showWarnings = FALSE)
dir.create(file.path(run$dir, "ldscores"), showWarnings = FALSE)
message("[06] run ", run$run_id, " opened")
cat(run$run_id, "\n", sep = "")
