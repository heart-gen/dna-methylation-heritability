#!/usr/bin/env Rscript
#### 08 Stage 02 -- tier 2: the one identified region difference ####
##
## Usage:
##   Rscript _h/02_identified_difference.R --run-id rdg-AA-crossregion-YYYYMMDD
##
## DLPFC vs hippocampus is the ONLY clean region contrast in this dataset. Both
## regions sit inside sequencing batches 1-2, so the comparison does not cross
## the batch boundary that makes every caudate contrast uninterpretable
## (AGENTS.md 8.1). That is the whole reason this tier exists and is restricted
## to two regions.
##
## Caudate is fitted, written and surfaced -- in its own columns, flagged, and
## excluded from the claim. That is the Module 04 mechanism
## (config/repeat_annotations.yml:interpretation.technically_confounded_regions)
## reused rather than reinvented: removing a region from a claim while keeping
## its estimate visible is a RESTRICTION of the evidence base, and dropping the
## rows instead would hide that a restriction was applied at all.
##
## Strict conjunction: a difference is claimed only when the primary contrast
## and every sensitivity agree. A difference that appears under one
## specification and not another is reported as a non-claim, not as a finding
## with caveats.

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
out_dir <- file.path(run_dir, "results")
tier <- "identified_difference"

contrast <- as.character(config_get(cfg, "identified_difference.contrast"))
if (length(contrast) != 2L) {
    stop("identified_difference.contrast must name exactly two regions")
}
confounded <- as.character(
    config_get(cfg, "interpretation.technically_confounded_regions.all_outcomes"))
alpha <- as.numeric(config_get(cfg, "identified_difference.alpha"))
fdr_method <- as.character(config_get(cfg, "identified_difference.fdr_method"))
strict <- isTRUE(config_get(cfg, "identified_difference.require_strict_conjunction"))

## The contrast must not include a region this config itself calls confounded.
## A config edit that put caudate into tier 2 would defeat the tier's premise,
## so refuse rather than compute.
if (length(intersect(contrast, confounded))) {
    stop("identified_difference.contrast includes a technically confounded ",
         "region (", paste(intersect(contrast, confounded), collapse = ", "),
         "). Tier 2 exists because it avoids the batch boundary.")
}

tests <- as.data.table(fread(file.path(out_dir, "cross-region-tests.tsv")))
if (!nrow(tests)) stop("Stage 01 produced no tests; run it first")

## --------------------------------------------------------------- the contrast
##
## A two-sample z that ADDS the per-region variances. This is the only place in
## the module where a DIFFERENCE of levels is computed, and it is legitimate here
## precisely because the two statistics come from the same analysis in the same
## units -- unlike the Module 02 score, which is a within-cell rank (AGENTS.md
## 7.6).
##
## The two regions are NOT disjoint donor sets: dlpfc and hippocampus share 115
## of 118 donors (Jaccard 0.96), so the estimates are positively correlated and
## Var(a - b) = Var(a) + Var(b) - 2Cov(a, b) is SMALLER than what this z uses.
## The test is therefore conservative -- it under-rejects -- which is the safe
## direction for a tier whose licence is "genuine regional heterogeneity". A
## paired donor bootstrap is the estimator that would use the covariance;
## 09b_aging_application/_h/05_cross_region_concordance.R implements one for the
## same pair. Replacing this z is a scope change, not a bug fix, and the earlier
## claim of independence in this comment was simply wrong.
## Keyed the same way Stage 01 keys: analysis_set distinguishes Module 04's
## primary fit from its four sensitivities, and omitting it would collapse five
## distinct tests into one cell of the cast.
##
## `outcome_role` travels with the test because tier 2 must count the way tier 1
## counts. It is a property of the test, not of the region, so it is carried on the
## LHS; a role that varied by region would split a test into two rows and silently
## double a difference, so the invariant is checked rather than assumed.
role_by_test <- unique(tests[, .(test_id, outcome_role)])
dupe_role <- role_by_test[, .N, by = test_id][N > 1L]
if (nrow(dupe_role)) {
    stop("outcome_role varies by region for ", nrow(dupe_role), " test(s), e.g. ",
         dupe_role$test_id[[1]], ". The role identifies the test, so a ",
         "region-dependent role would split one test into two contrast rows.")
}
wide <- dcast(tests[region %in% c(contrast, confounded)],
              analysis + analysis_set + outcome_role + outcome + predictor +
                  test_id ~ region,
              value.var = c("estimate", "se", "p", "q", "n"))

a <- contrast[[1]]; b <- contrast[[2]]
ea <- paste0("estimate_", a); eb <- paste0("estimate_", b)
sa <- paste0("se_", a);       sb <- paste0("se_", b)
for (col in c(ea, eb, sa, sb)) {
    if (!col %in% names(wide)) {
        stop("No ", col, " column; tier 2 needs both contrast regions in ",
             "Stage 01's output")
    }
}

wide[, delta := get(ea) - get(eb)]
wide[, delta_se := sqrt(get(sa)^2 + get(sb)^2)]
wide[, delta_z := fifelse(is.finite(delta_se) & delta_se > 0,
                          delta / delta_se, NA_real_)]
wide[, delta_p := 2 * stats::pnorm(-abs(delta_z))]
wide[, contrast := paste(a, "minus", b)]

## Only complete pairs are testable. Say so in a column rather than dropping
## them, so the denominator stays auditable.
wide[, testable := is.finite(delta) & is.finite(delta_z)]
wide[testable == TRUE, delta_q := stats::p.adjust(delta_p, method = fdr_method)]

## ------------------------------------------------------- strict conjunction
##
## Sensitivities: the difference must survive each one, not just the primary
## fit. `both_nominal` guards the degenerate case where a "difference" is
## really one region estimating something and the other estimating nothing --
## that is a power contrast, not regional heterogeneity.
##
## BOTH regions must be nominal, so the operator is `&`. It was `|` until
## 2026-09-23, which made the check vacuous: a pair in which only one region
## reaches p < 0.05 is exactly the degenerate case the paragraph above says this
## guard exists to exclude, and an OR admitted it. The bug was not academic --
## it was the sole reason tier 2 reported a claimed difference at all
## (expression_abc x any_meqtl_support: hippocampus near-complete separation on
## 8 coupled VMRs of 244, DLPFC p = 0.0857). Under `&` the count is zero of 154
## testable pairs, which is the honest reading and is what
## writing-notes/WGBS_BATCH_REGION_CONFOUNDING.md and 07/README.md:107 already
## say about that outcome independently.
wide[, sens_fdr_significant := is.finite(delta_q) & delta_q < alpha]
wide[, sens_both_nominal := is.finite(get(paste0("p_", a))) &
         is.finite(get(paste0("p_", b))) &
         (get(paste0("p_", a)) < 0.05 & get(paste0("p_", b)) < 0.05)]
wide[, sens_direction_opposed := sign(get(ea)) != sign(get(eb))]
wide[, sens_magnitude := abs(delta) > pmax(get(sa), get(sb))]

sens_cols <- c("sens_fdr_significant", "sens_both_nominal", "sens_magnitude")
wide[, n_sensitivities_passed := rowSums(
    as.matrix(.SD == TRUE), na.rm = TRUE), .SDcols = sens_cols]
wide[, difference_claimed := if (strict) {
    testable & n_sensitivities_passed == length(sens_cols)
} else {
    testable & sens_fdr_significant
}]
## ------------------------------------------------ what counts as a difference
##
## `difference_claimed` above answers a STATISTICAL question about one row: did
## this contrast survive the conjunction. Whether that row is a countable finding
## is a separate question, and tier 2 was not asking it. Stage 01 asks it for tier
## 1 and its comments state why: Module 04's analysis_sets are the same tests
## refit, so counting all of them "would report one test five or six times and
## inflate the primary deliverable roughly fivefold", and "a negative control is
## not support". Tier 2 counted every analysis_set and every outcome_role, which is
## the inflation stage 01 refuses by design -- and the config points the same way,
## naming Module 04 as `template_module` for its "primary claim on the identified
## contrast".
##
## So the headline is the primary arm of the prespecified claim family, exactly as
## in stage 01. The unfiltered count stays in the summary as the auditable
## denominator, and every row stays in the table with its own verdict; nothing is
## dropped, only counted correctly.
wide[, is_primary := analysis_set == "primary"]
wide[, in_claim_family := outcome_role %in% c("bh_family", "prespecified_family")]
wide[, is_negative_control := outcome_role == "specificity_control"]
wide[, difference_claimed_primary_claim_family :=
         difference_claimed & is_primary & in_claim_family]
## A control differing between regions is not a regional-heterogeneity finding; it
## is a note about the control. Counted separately so it can never be added in.
wide[, difference_claimed_negative_control :=
         difference_claimed & is_primary & is_negative_control]
wide[, tier := tier]
wide[, licenses := trimws(config_get(cfg, paste0("tiers.", tier, ".licenses")))]

write_atomic(wide[order(delta_p)],
             file.path(out_dir, "identified-difference.tsv"))

## ------------------------------------------- tier 4: caudate, set aside
##
## Fitted, written, surfaced, and excluded from every claim. The estimate is
## kept so that a reader who asks "what about caudate" gets a number and a
## reason rather than silence.
descriptive <- tests[region %in% confounded]
descriptive[, tier := "descriptive_only"]
descriptive[, licenses := trimws(
    config_get(cfg, "tiers.descriptive_only.licenses"))]
descriptive[, set_aside_reason :=
                "region perfectly confounded with sequencing batch (AGENTS.md 8.1)"]
descriptive[, supports_any_claim := FALSE]
write_atomic(descriptive[order(analysis, analysis_set, outcome, predictor)],
             file.path(out_dir, "descriptive-confounded-regions.tsv"))

summary_dt <- data.table(
    tier = tier,
    contrast = paste(a, "minus", b),
    batches = paste(config_get(cfg, "identified_difference.batches"),
                    collapse = ","),
    strict_conjunction = strict,
    n_pairs = nrow(wide),
    n_testable = wide[testable == TRUE, .N],
    n_fdr_significant = wide[sens_fdr_significant == TRUE, .N],
    ## Every row that survived the statistical conjunction, all analysis_sets and
    ## all roles. Reported so the headline's denominator is auditable, never as the
    ## headline itself.
    n_difference_claimed = wide[difference_claimed == TRUE, .N],
    n_difference_claimed_primary = wide[difference_claimed == TRUE &
                                            is_primary == TRUE, .N],
    ## THE HEADLINE: primary arm, prespecified claim family. One finding counted
    ## once, and a negative control never counted as heterogeneity.
    n_difference_claimed_primary_claim_family =
        wide[difference_claimed_primary_claim_family == TRUE, .N],
    n_difference_claimed_negative_control =
        wide[difference_claimed_negative_control == TRUE, .N],
    n_claim_family_pairs = wide[is_primary == TRUE & in_claim_family == TRUE, .N],
    alpha = alpha,
    fdr_method = fdr_method,
    confounded_regions_set_aside = paste(confounded, collapse = ","),
    n_descriptive_rows_retained = nrow(descriptive))
write_atomic(summary_dt, file.path(out_dir, "identified-difference-summary.tsv"))

print(summary_dt)
message("[tier2] ", summary_dt$n_difference_claimed_primary_claim_family, " of ",
        summary_dt$n_claim_family_pairs, " primary claim-family ", a, "-vs-", b,
        " differences survive strict conjunction")
message("  Denominators: ", summary_dt$n_testable, " testable pairs across all ",
        "analysis_sets and roles, of which ", summary_dt$n_difference_claimed,
        " survive the conjunction and ", summary_dt$n_difference_claimed_primary,
        " are in the primary arm. The headline counts the primary arm of the ",
        "prespecified claim family only: a sensitivity refit is not a second ",
        "difference, a rescaled outcome is not a second finding, and a ",
        "specificity control differing between regions (",
        summary_dt$n_difference_claimed_negative_control,
        ") is a note about the control, not regional heterogeneity.")
message("  ", nrow(descriptive), " ", paste(confounded, collapse = "/"),
        " rows retained as DESCRIPTIVE ONLY: fitted and surfaced, excluded ",
        "from every claim.")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
