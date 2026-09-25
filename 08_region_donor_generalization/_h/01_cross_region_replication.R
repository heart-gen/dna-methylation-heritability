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
if (run_is_sealed(manifest)) {
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
## `analysis_set` is part of the KEY, not decoration. Module 04 publishes its
## primary fit plus four sensitivities in every region (primary,
## high_mappability, exclude_segdups, low_cell_composition,
## adjust_cell_composition) -- 5 sets x 28 outcome/predictor pairs = 140 rows --
## and, in caudate only, the conditional adjust_cell_composition_scmd arm, which
## is why the set count is read from the table rather than assumed. Keying on
## outcome+predictor alone makes every Module 04 test appear 15 times instead
## of 3, so the "complete in all three regions" filter drops the entire module
## -- the primary biological analysis -- without a word. Modules 05 and 07
## publish a single set and get "primary".
##
## `outcome_role` is carried for the same reason. Module 04 publishes five
## roles, and they are not interchangeable evidence:
##   bh_family            the prespecified FDR-controlled family -- the claim
##   specificity_control  SHOULD be null; replicating is a WARNING, not support
##   secondary_scale      the same outcome on another scale -- double-counting
##   complementary_contrast / descriptive  context, not claim
## A headline that pooled all five would count one finding several times and
## would count a failed negative control as a success. Modules 05 and 07 have
## no role column and are treated as their own prespecified family.
COMMON <- c("region", "analysis", "analysis_set", "outcome_role", "outcome",
            "predictor", "estimate", "se", "p", "q", "n")

harvest <- function(module, region, file, mapping, analysis_label) {
    f <- run_results(module, region, file)
    if (!file.exists(f)) {
        stop("Upstream result missing: ", f,
             "\n  Tier 1 assembles accepted results; it cannot recompute them.")
    }
    dt <- as.data.table(fread(f))
    ## A row an upstream marked NOT fitted is not a test. Module 04 records its
    ## conditional scMD composition arm in every region -- with estimates only
    ## where the integration gate passes (AGENTS.md 7.4) -- so that the absence is
    ## visible rather than silent. Harvesting those rows would make the arm look
    ## complete in all three regions and then fail direction consistency on an NA
    ## estimate, i.e. it would report a sensitivity that was never run as a
    ## sensitivity that did not replicate.
    if ("arm_fitted" %in% names(dt)) {
        n_unfitted <- sum(!(dt$arm_fitted %in% TRUE))
        if (n_unfitted > 0L) {
            message("[08] ", module, "/", region, ": dropped ", n_unfitted,
                    " row(s) marked arm_fitted = FALSE")
        }
        dt <- dt[dt$arm_fitted %in% TRUE]
    }
    ## Every mapping entry names a COLUMN except `nuisance_terms`, which is a
    ## regex matched against the predictor column's VALUES. Including it here
    ## makes the check look for a column literally named
    ## "^(\\(Intercept\\)|cpg_density|...)$" and refuse the Module 05 table.
    cols <- mapping[setdiff(names(mapping), "nuisance_terms")]
    missing <- setdiff(unlist(cols), names(dt))
    if (length(missing)) {
        stop(module, " ", basename(f), " lacks: ",
             paste(missing, collapse = ", "))
    }
    out <- data.table(
        region    = region,
        analysis  = analysis_label,
        analysis_set = if (!is.null(mapping$analysis_set))
                           as.character(dt[[mapping$analysis_set]])
                       else "primary",
        outcome_role = if (!is.null(mapping$outcome_role)) {
                           as.character(dt[[mapping$outcome_role]])
                       } else if (!is.null(mapping$nuisance_terms)) {
                           ## A model's covariate and intercept rows are not
                           ## scientific tests. Module 05 publishes its whole
                           ## design matrix, so without this the tier-1
                           ## denominator would include (Intercept).
                           ifelse(grepl(mapping$nuisance_terms,
                                        as.character(dt[[mapping$predictor]])),
                                  "nuisance", "prespecified_family")
                       } else "prespecified_family",
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
         mapping = list(analysis_set = "analysis_set",
                        outcome_role = "outcome_role", outcome = "outcome",
                        predictor = "predictor",
                        estimate = "estimate", se = "se", p = "p", q = "q",
                        n = "n")),
    list(module = "05_cpg_meqtl_burden",
         file = "burden-primary-model.tsv", analysis = "meqtl_burden",
         ## The claim is the predictability gradient. The intercept, CpG
         ## density, VMR length and mean methylation are the covariates it is
         ## adjusted for, not findings, so they are marked nuisance and stay
         ## out of the headline denominator while remaining in the table.
         mapping = list(outcome = "model", predictor = "term",
                        nuisance_terms = "^(\\(Intercept\\)|cpg_density|log\\(vmr_length\\)|mean_methylation)$",
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

## A test is identified by what it tests, not by where it ran. A duplicated key
## WITHIN one region means the key is incomplete for that source -- the exact
## defect analysis_set fixes -- so stop rather than let it inflate n_regions.
tests[, test_id := paste(analysis, analysis_set, outcome, predictor,
                         sep = "::")]
dupe <- tests[, .N, by = .(region, test_id)][N > 1L]
if (nrow(dupe)) {
    stop("Duplicated test keys within a region (", nrow(dupe), " cases, e.g. ",
         dupe$test_id[[1]], " x", dupe$N[[1]], "). The harvest key is ",
         "incomplete for that source; add the distinguishing column to its ",
         "mapping rather than letting it inflate n_regions.")
}
tests[, direction := fifelse(estimate > 0, "up",
                             fifelse(estimate < 0, "down", "flat"))]
tests[, tier := tier]
write_atomic(tests[order(analysis, analysis_set, outcome, predictor, region)],
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
), by = .(analysis, analysis_set, outcome_role, outcome, predictor,
             test_id)]

## Hoisted OUT of the data.table `[`: per_test carries its own `regions` column
## (the comma-joined region list), so the bare name inside `[` resolves to that
## COLUMN and length(regions) is the ROW COUNT, not 3. The test silently became
## `n_regions == 154`, false for every row, and tier 1 reported zero tests
## complete and zero replications while every test was in fact present in all
## three regions. Same defect class as the `cohort`/`region` shadowing in
## gates.R:require_accepted_upstream() and the `cell` shadowing in stage 04.
n_regions_expected <- length(regions)
stopifnot(n_regions_expected >= 1L, !is.na(n_regions_expected))
per_test[, complete_across_regions := n_regions == n_regions_expected]
per_test[, direction_consistent := pmax(n_up, n_down) == n_regions]
## The replication call. Consistent direction AND nominal support in at least
## `min_regions` regions: direction alone would let three null estimates that
## happen to share a sign count as a replication.
per_test[, replicated := complete_across_regions & direction_consistent &
             n_nominal >= min_regions]

## -------------------------------------------- primary vs sensitivity sets
##
## Module 04's analysis_sets are the SAME 28 tests refit under its sensitivities,
## so counting all of them as replications would report one test five or six times
## and inflate the primary deliverable roughly fivefold. The
## headline count is therefore PRIMARY ONLY, and the sensitivities are used the
## way Module 04 itself uses them: as a strict conjunction the primary claim
## must survive, not as extra evidence.
per_test[, is_primary := analysis_set == "primary"]
sens <- per_test[is_primary == FALSE,
                 .(n_sensitivity_sets = .N,
                   n_sensitivity_replicated = sum(replicated)),
                 by = .(analysis, outcome_role, outcome, predictor)]
per_test <- merge(per_test, sens,
                  by = c("analysis", "outcome_role", "outcome", "predictor"),
                  all.x = TRUE)
per_test[is.na(n_sensitivity_sets), `:=`(n_sensitivity_sets = 0L,
                                         n_sensitivity_replicated = 0L)]
## A primary test replicates STRICTLY when it replicates and so does every
## sensitivity refit of it. A test with no sensitivities passes trivially, and
## the n_sensitivity_sets column is what tells a reader which case they have.
per_test[, replicated_strict := replicated & is_primary &
             n_sensitivity_replicated == n_sensitivity_sets]
## The claim family: the prespecified, FDR-controlled tests only.
per_test[, in_claim_family := outcome_role %in% c("bh_family",
                                                  "prespecified_family")]
per_test[, is_negative_control := outcome_role == "specificity_control"]

## A specificity control is NOT judged by whether it replicates -- it is judged
## by WHICH WAY. Replicating in the opposite direction to the claim family is
## the control doing its job: it says the association is specific to the tracks
## claimed rather than a property of any annotation. Replicating in the SAME
## direction is the finding that undercuts the tier. Scoring both as "a control
## replicated" would flag the healthy case and read the two identically.
claim_sign <- per_test[in_claim_family == TRUE & is_primary == TRUE,
                       sign(sum(sign(n_up - n_down)))]
per_test[, claim_family_direction := fifelse(claim_sign > 0, "up",
                                             fifelse(claim_sign < 0, "down",
                                                     "mixed"))]
per_test[, control_direction := fifelse(n_up > n_down, "up",
                                        fifelse(n_down > n_up, "down", "mixed"))]
per_test[, control_opposes_claim := is_negative_control &
             claim_family_direction != "mixed" &
             control_direction != "mixed" &
             control_direction != claim_family_direction]
per_test[, control_tracks_claim := is_negative_control & replicated &
             control_direction == claim_family_direction]
per_test[, tier := tier]
per_test[, licenses := trimws(config_get(cfg, paste0("tiers.", tier, ".licenses")))]
write_atomic(per_test[order(-in_claim_family, -replicated, analysis,
                              analysis_set, outcome, predictor)],
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
    n_rows_harvested = nrow(tests),
    ## Distinct tests, all analysis_sets. Reported so the denominator of the
    ## primary counts below is auditable, never as the headline.
    n_tests_all_sets = nrow(per_test),
    n_tests_complete_all_sets = per_test[complete_across_regions == TRUE, .N],
    ## THE HEADLINE. Primary analysis_set, prespecified claim family only, so
    ## one finding is counted once and a negative control is never counted as
    ## support.
    n_claim_tests = per_test[is_primary == TRUE & in_claim_family == TRUE, .N],
    n_claim_complete = per_test[is_primary == TRUE & in_claim_family == TRUE &
                                    complete_across_regions == TRUE, .N],
    n_claim_direction_consistent = per_test[is_primary == TRUE &
                                                in_claim_family == TRUE &
                                                complete_across_regions == TRUE &
                                                direction_consistent == TRUE, .N],
    n_claim_replicated = per_test[is_primary == TRUE & in_claim_family == TRUE &
                                      replicated == TRUE, .N],
    n_claim_replicated_strict = per_test[replicated_strict == TRUE &
                                             in_claim_family == TRUE, .N],
    ## Reported alongside, never added in.
    n_secondary_scale_replicated = per_test[is_primary == TRUE &
                                                outcome_role == "secondary_scale" &
                                                replicated == TRUE, .N],
    ## Controls, split by direction. `tracks_claim` is the red flag;
    ## `opposes_claim` is the control working as designed.
    claim_family_direction = per_test$claim_family_direction[[1]],
    n_negative_controls = per_test[is_primary == TRUE &
                                       is_negative_control == TRUE, .N],
    n_controls_opposing_claim = per_test[is_primary == TRUE &
                                             control_opposes_claim == TRUE, .N],
    n_controls_tracking_claim = per_test[is_primary == TRUE &
                                             control_tracks_claim == TRUE, .N],
    n_nuisance_terms_excluded = per_test[is_primary == TRUE &
                                             outcome_role == "nuisance", .N],
    n_tests_primary = per_test[is_primary == TRUE, .N],
    n_sensitivity_sets_max = max(per_test$n_sensitivity_sets),
    min_regions_for_replication = min_regions,
    ## Carried so no consumer has to know the rule to read the table.
    raw_score_comparison_emitted = FALSE,
    agreement_basis = config_get(cfg, "cross_region_replication.agreement_basis"))
write_atomic(summary_dt, file.path(out_dir, "cross-region-summary.tsv"))

print(summary_dt)
message("[tier1] claim family: ", summary_dt$n_claim_replicated, " of ",
        summary_dt$n_claim_complete,
        " prespecified tests complete in all three regions replicate; ",
        summary_dt$n_claim_replicated_strict,
        " also survive every sensitivity refit")
message("  Denominators: ", summary_dt$n_tests_all_sets,
        " tests across all analysis_sets, ", summary_dt$n_tests_primary,
        " in the primary set, ", summary_dt$n_claim_tests,
        " in the prespecified claim family. The headline counts the claim ",
        "family only: a sensitivity refit is not a second test, a rescaled ",
        "outcome is not a second finding, and a negative control is not ",
        "support.")
if (summary_dt$n_controls_tracking_claim > 0L) {
    message("  WARNING: ", summary_dt$n_controls_tracking_claim, " of ",
            summary_dt$n_negative_controls, " specificity controls replicate ",
            "in the SAME direction as the claim family (",
            summary_dt$claim_family_direction, "). That undercuts the tier: ",
            "the association is not specific to the tracks claimed.")
}
if (summary_dt$n_controls_opposing_claim > 0L) {
    message("  ", summary_dt$n_controls_opposing_claim, " of ",
            summary_dt$n_negative_controls, " specificity controls run ",
            "OPPOSITE to the claim family. That is the control working as ",
            "designed -- evidence of specificity, not a failure.")
}
message("  Robustness across technical AND regional contexts. Region is ",
        "confounded with sequencing batch (AGENTS.md 8.1): a successful ",
        "replication is more compelling for crossing it, a difference is not ",
        "more interpretable.")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
