#!/usr/bin/env Rscript
#### 10_environmental_exploratory -- open a run ####
##
## Usage:
##   Rscript _h/00_new_run.R --cohort AA --region caudate [--allow-unlocked]
##
## Three upstreams are consumed: 01 for the corrected VMR boundaries and the
## per-locus methylation phenotypes, 02 for the relative local-control score
## that replaces v1's retired h2_category, and 04 for the technical covariates
## the axis models adjust on. All must describe the same vmr_set_id, or the
## module would join an exposure association, a rank and a GC content computed
## for different loci sharing an ID.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "10_environmental_exploratory"
MODULE_TAG <- "env"

opts <- parse_v2_args(require = c("cohort", "region"))
allow_unlocked <- isTRUE(opts$allow_unlocked)

env <- load_config("environmental")
assert_locked(list(environmental = env), allow_unlocked = allow_unlocked)

## AGENTS.md 3 and 2.3. v1's endpoints were all read against h2_category, built
## from h2_unscaled and r_squared_cv. Checking the predictor here means the ban
## is enforced before a run directory exists, not after a SLURM array.
if (!identical(env$testing$architecture_predictor,
               "local_snp_contribution_score_z")) {
    stop("config/environmental.yml must use local_snp_contribution_score_z as ",
         "the architecture predictor. Legacy predictability and absolute PVE ",
         "are banned (AGENTS.md 3), and the heritable/non-heritable grouping ",
         "v1 used is banned outright (AGENTS.md 2.3).")
}
if (!isTRUE(env$testing$separate_fdr_family)) {
    stop("testing.separate_fdr_family must be true: one BH family per ",
         "exposure, never pooled (AGENTS.md 10.3).")
}
if (!isTRUE(env$interpretation$exploratory_supplement_only)) {
    stop("interpretation.exploratory_supplement_only must be true. AGENTS.md ",
         "2.3 permits exposure results only as descriptive or sensitivity ",
         "analyses in the supplement.")
}
if (isTRUE(env$interpretation$cross_region_comparison_allowed)) {
    stop("interpretation.cross_region_comparison_allowed must be false: ",
         "region is perfectly confounded with sequencing batch (AGENTS.md 8.1).")
}
for (key in c("min_minor_class_n", "max_missing_frac")) {
    if (is.null(env$eligibility[[key]])) {
        stop("eligibility.", key, " is unset. The gate is the concrete form of ",
             "the power limitation that demoted this analysis; it must be ",
             "prespecified, not chosen after the counts are seen.")
    }
}

upstreams <- list(
    vmr_catalog            = "01_vmr_catalog",
    local_genetic_variance = "02_local_genetic_variance",
    repeat_architecture    = "04_repeat_repressive_architecture"
)
accepted <- lapply(upstreams, require_accepted_upstream,
                   cohort = opts$cohort, region = opts$region,
                   allow_unaccepted = allow_unlocked)

sets <- unlist(lapply(accepted, function(a) a$vmr_set_id))
sets <- unique(sets[!is.na(sets)])
if (length(sets) > 1) {
    stop("vmr_set_id mismatch across upstreams: ", paste(sets, collapse = " vs "),
         "\n  01, 02 and 04 must describe the same VMR set (AGENTS.md 6).")
}

if (!is.null(opts$run_id) && !allow_unlocked) {
    stop("--run-id may only be given for a smoke run (--allow-unlocked); ",
         "production run IDs are derived, not chosen.")
}

eligible_exposures <- c(names(env$exposures$binary),
                        names(env$exposures$categorical))

run <- new_run(
    module = MODULE_TAG, cohort = opts$cohort, region = opts$region,
    module_root = file.path(repo_root(), MODULE),
    run_id = opts$run_id,
    vmr_set_id = if (length(sets)) sets[1] else NA_character_,
    upstream = stats::setNames(
        lapply(accepted, function(a) a$run_id %||% NA_character_),
        paste0(names(accepted), "_run_id")),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        config_environmental_sha256 = attr(env, "config_sha256"),
        catalog_cohort = opts$catalog_cohort,
        estimation_group = opts$estimation_group,
        candidate_exposures = paste(eligible_exposures, collapse = ","),
        min_minor_class_n = as.character(env$eligibility$min_minor_class_n),
        max_missing_frac = as.character(env$eligibility$max_missing_frac),
        fdr_alpha = as.character(env$testing$fdr_alpha),
        primary_axis_model = env$testing$primary_axis_model,
        exploratory_supplement_only = "TRUE"
    )
)

for (d in c("results", "results/per-chrom", "logs")) {
    dir.create(file.path(run$dir, d), showWarnings = FALSE, recursive = TRUE)
}
message("[10] run ", run$run_id, " opened (candidate exposures: ",
        paste(eligible_exposures, collapse = ", "), ")")
cat(run$run_id, "\n", sep = "")
