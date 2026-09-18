#!/usr/bin/env Rscript
#### 08 Stage 01 -- tier 1: cross-region replication ####
##
## Usage:
##   Rscript _h/01_cross_region_replication.R --run-id rdg-AA-crossregion-YYYYMMDD
##
## The PRIMARY deliverable (AGENTS.md 7.7 tier 1). It assembles results that are
## already accepted per region and asks whether they agree, region by region. It
## refits nothing.
##
## What "agreement" may mean here is tightly bounded. AGENTS.md 7.6 forbids
## comparing the Module 02 score across regions at all, so agreement is measured
## on the DIRECTION and the RANK of each already-published test statistic, never
## on its level. Two regions agreeing that a predictor acts in the same
## direction with the same sign of effect is a replication; their coefficients
## being numerically close is not something this module is allowed to assert.
##
## Read the result as robustness across technical AND regional contexts,
## because region is perfectly confounded with sequencing batch (AGENTS.md 8.1).
## A successful replication is more compelling for crossing that boundary. A
## failure is NOT thereby more interpretable, and the config's
## tiers.cross_region_replication.forbid_phrasings names the sentences that get
## this backwards.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "08_region_donor_generalization"

opts <- parse_v2_args(require = c("run_id"))
run_dir <- file.path(V2_ROOT, MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)
manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]
    if (length(v) != 1L) stop("Run manifest lacks unique field: ", f)
    as.character(v[[1L]])
}
if (nzchar(manifest$value[manifest$field == "finished_at"][1] %||% "")) {
    stop("Run is sealed and immutable: ", opts$run_id)
}

cfg <- load_config("region_donor_generalization")
cohort <- mval("cohort")
regions <- strsplit(mval("regions"), ",", fixed = TRUE)[[1]]
tier <- "cross_region_replication"
min_regions <- as.integer(config_get(cfg, "cross_region_replication.min_regions_for_replication"))
out_dir <- file.path(run_dir, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

upstream_run <- function(module, region) {
    mval(paste0("upstream_", gsub("[|.]", "_", paste(module, region, sep = "|"))))
}
run_results <- function(module, region, file) {
    file.path(V2_ROOT, module, "_m", "runs", upstream_run(module, region),
              "results", file)
}

## --------------------------------------------------------------- harvesting
##
## Each upstream publishes a long table of tests with estimate/se/p/q. We take
## the columns that are common to all of them and nothing else, so a schema
## change in one module cannot silently reshape the replication table.
COMMON <- c("region", "analysis", "outcome", "predictor", "estimate", "se",
            "p", "q", "n")

harvest <- function(module, region, file, mapping, analysis_label) {
    f <- run_results(module, region, file)
    if (!file.exists(f)) {
        stop("Upstream result missing: ", f,
             "\n  Tier 1 assembles accepted results; it cannot recompute them.")
    }
    dt <- as.data.table(fread(f))
    missing <- setdiff(unlist(mapping), names(dt))
    if (length(missing)) {
        stop(module, " ", basename(f), " lacks: ",
             paste(missing, collapse = ", "))
    }
    out <- data.table(
        region    = region,
        analysis  = analysis_label,
        outcome   = as.character(dt[[mapping$outcome]]),
        predictor = as.character(dt[[mapping$predictor]]),
        estimate  = suppressWarnings(as.numeric(dt[[mapping$estimate]])),
        se        = suppressWarnings(as.numeric(dt[[mapping$se]])),
        p         = suppressWarnings(as.numeric(dt[[mapping$p]])),
        q         = if (!is.null(mapping$q))
                        suppressWarnings(as.numeric(dt[[mapping$q]]))
                    else NA_real_,
        n         = if (!is.null(mapping$n))
                        suppressWarnings(as.numeric(dt[[mapping$n]]))
                    else NA_real_,
        source_module = module,
        source_run    = upstream_run(module, region))
    out[is.finite(estimate)]
}

specs <- list(
    list(module = "04_repeat_repressive_architecture",
         file = "association-results.tsv", analysis = "repeat_architecture",
         mapping = list(outcome = "outcome", predictor = "predictor",
                        estimate = "estimate", se = "se", p = "p", q = "q",
                        n = "n")),
    list(module = "05_cpg_meqtl_burden",
         file = "burden-primary-model.tsv", analysis = "meqtl_burden",
         mapping = list(outcome = "model", predictor = "term",
                        estimate = "estimate", se = "se", p = "p",
                        n = "n_vmrs")),
    list(module = "07_transcription_splicing_coupling",
         file = "coupling-tests.tsv", analysis = "expression_coupling",
         mapping = list(outcome = "modality", predictor = "predictor",
                        estimate = "estimate", se = "se", p = "p", q = "q",
                        n = "n"))
)

tests <- rbindlist(lapply(regions, function(re) {
    rbindlist(lapply(specs, function(s)
        harvest(s$module, re, s$file, s$mapping, s$analysis)), use.names = TRUE)
}), use.names = TRUE)

## A test is identified by what it tests, not by where it ran.
tests[, test_id := paste(analysis, outcome, predictor, sep = "::")]
tests[, direction := fifelse(estimate > 0, "up",
                             fifelse(estimate < 0, "down", "flat"))]
tests[, tier := tier]
write_atomic(tests[order(analysis, outcome, predictor, region)],
             file.path(out_dir, "cross-region-tests.tsv"))

## ------------------------------------------------------------- replication
##
## Direction agreement, region by region, on tests present in every region.
## Restricting to the complete set is not cosmetic: a test present in two
## regions and absent in the third would otherwise be scored as replicating
## everywhere it was measured, which reads as stronger than it is.
per_test <- tests[, .(
    n_regions      = .N,
    regions        = paste(sort(region), collapse = ","),
    n_up           = sum(direction == "up"),
    n_down         = sum(direction == "down"),
    n_nominal      = sum(is.finite(p) & p < 0.05),
    n_fdr          = sum(is.finite(q) & q < 0.05),
    min_p          = if (any(is.finite(p))) min(p, na.rm = TRUE) else NA_real_,
    estimate_range = diff(range(estimate))
), by = .(analysis, outcome, predictor, test_id)]

per_test[, complete_across_regions := n_regions == length(regions)]
per_test[, direction_consistent := pmax(n_up, n_down) == n_regions]
## The replication call. Consistent direction AND nominal support in at least
## `min_regions` regions: direction alone would let three null estimates that
## happen to share a sign count as a replication.
per_test[, replicated := complete_across_regions & direction_consistent &
             n_nominal >= min_regions]
per_test[, tier := tier]
per_test[, licenses := trimws(config_get(cfg, paste0("tiers.", tier, ".licenses")))]
write_atomic(per_test[order(-replicated, analysis, outcome, predictor)],
             file.path(out_dir, "cross-region-replication.tsv"))

## ------------------------------------------------- rank agreement per region
##
## Ordering agreement across regions, computed WITHIN an analysis so that
## incommensurable statistics are never ranked against each other. This is the
## rank half of `agreement_basis: rank_and_direction`; the level is never used.
pairs <- combn(regions, 2, simplify = FALSE)
rank_rows <- rbindlist(lapply(unique(tests$analysis), function(an) {
    rbindlist(lapply(pairs, function(pr) {
        a <- tests[analysis == an & region == pr[[1]]]
        b <- tests[analysis == an & region == pr[[2]]]
        shared <- intersect(a$test_id, b$test_id)
        if (length(shared) < 3L) {
            return(data.table(analysis = an, region_a = pr[[1]],
                              region_b = pr[[2]], n_shared = length(shared),
                              spearman_estimate = NA_real_,
                              direction_agreement = NA_real_))
        }
        av <- a[match(shared, test_id)]
        bv <- b[match(shared, test_id)]
        data.table(
            analysis = an, region_a = pr[[1]], region_b = pr[[2]],
            n_shared = length(shared),
            spearman_estimate = stats::cor(av$estimate, bv$estimate,
                                           method = "spearman"),
            direction_agreement = mean(av$direction == bv$direction))
    }), use.names = TRUE)
}), use.names = TRUE)
rank_rows[, tier := tier]
write_atomic(rank_rows, file.path(out_dir, "cross-region-rank-agreement.tsv"))

## ----------------------------------------------------------------- summary
summary_dt <- data.table(
    tier = tier,
    cohort = cohort,
    regions = paste(regions, collapse = ","),
    n_tests_total = nrow(tests),
    n_tests_complete = per_test[complete_across_regions == TRUE, .N],
    n_direction_consistent = per_test[complete_across_regions == TRUE &
                                          direction_consistent == TRUE, .N],
    n_replicated = per_test[replicated == TRUE, .N],
    min_regions_for_replication = min_regions,
    ## Carried so no consumer has to know the rule to read the table.
    raw_score_comparison_emitted = FALSE,
    agreement_basis = config_get(cfg, "cross_region_replication.agreement_basis"))
write_atomic(summary_dt, file.path(out_dir, "cross-region-summary.tsv"))

print(summary_dt)
message("[tier1] ", summary_dt$n_replicated, " of ",
        summary_dt$n_tests_complete,
        " tests complete in all three regions replicate")
message("  Robustness across technical AND regional contexts. Region is ",
        "confounded with sequencing batch (AGENTS.md 8.1): a successful ",
        "replication is more compelling for crossing it, a difference is not ",
        "more interpretable.")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
