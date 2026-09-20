#!/usr/bin/env Rscript
#### 09 stage 15 -- cross-region axis concordance, and decision 2 ####
##
## Usage:
##   Rscript _h/15_cross_region_axis_concordance.R --cohort AA
##   Rscript _h/15_cross_region_axis_concordance.R --cohort AA --allow-unlocked
##
## PI 2026-09-18. Module 09's framing question is "does regional variation in
## local genetic control of methylation intersect schizophrenia-relevant
## regulatory biology, and is that relationship shared or region-dependent?"
## Every other stage of this module runs inside ONE region, so no stage could
## answer the "shared or region-dependent" half. This one does, and it is the
## only place decision 2 is resolved.
##
## H1, as locked (config independent_decisions.scz_application_retention):
##   SCZ-linked VMRs show LOWER local genetic-control scores, with statistical
##   support in AT LEAST TWO regions including AT LEAST ONE non-caudate region.
## A negative-but-nonsignificant region does NOT falsify H1.
##
## BIOLOGICAL UNIFORMITY IS NOT REQUIRED (PI 2026-09-18). Two questions are asked
## separately, and both are answered:
##
##   1. REGION-GENERAL ASSOCIATION -- is there a disorder-related association
##      that holds across regions? Decided by the >=2-regions rule above.
##   2. REGION-SPECIFIC DEPARTURE -- does any region genuinely depart from it?
##      Tested, not assumed away. A departure QUALIFIES the scope of claim 1; it
##      does not veto it (departure_vetoes_region_general_claim: false).
##
## Requiring the same direction everywhere would have done two wrong things at
## once: let one deviating region falsify a real region-general association, and
## bury a genuine region-specific finding by scoring it as failure.
##
## WHY CAUDATE CANNOT BE THE ONLY SUPPORT: it is perfectly confounded with
## sequencing batch (AGENTS.md 8.1), so a pattern resting on caudate alone is not
## separable from batch. It may contribute to concordance; it may not carry it.
##
## WHAT THIS STAGE MUST NOT DO: pool the three regions into one p-value. The
## regions share donors heavily -- 100 donors are in all three, DLPFC and
## hippocampus overlap in 115 of 118 (Jaccard 0.96), and the union is 168 donors
## -- so the three axis tests are NOT independent replicates and an
## inverse-variance meta-analysis would overstate precision by assuming they
## are. The endpoint is CONCORDANCE of sign and per-region significance. Donor
## overlap is computed and written out so the figure legend and the manuscript
## cannot describe this as independent replication.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09_schizophrenia_risk_application", "_h")),
    "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "09_schizophrenia_risk_application"
opts <- parse_v2_args(require = "cohort")
## --allow-unlocked is the repo's only valueless flag (00_shared/config.R:137);
## reuse it rather than inventing a second spelling for the same idea.
allow_unaccepted <- isTRUE(opts$allow_unlocked)

scz <- load_config("schizophrenia")
rule <- config_get(scz, "independent_decisions.scz_application_retention")
predictor <- scz$testing$architecture_predictor
alpha <- as.numeric(scz$testing$fdr_alpha)

dir_required   <- as.character(rule$axis_direction_required)
need_all_dir   <- isTRUE(rule$axis_direction_consistent_in_all_regions)
min_support    <- as.integer(rule$axis_min_regions_with_statistical_support)
support_outside <- as.character(rule$axis_requires_support_outside)

module_root <- file.path(repo_root(), MODULE)
out_dir <- file.path(module_root, "_m", "combined")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

regions <- as.character(config_get(load_config("cohorts"), "regions"))
if (length(regions) < 2L) stop("Need at least two regions for a concordance call")

## ------------------------------------------------------------- the three runs
##
## Acceptance is required by default (AGENTS.md 6): a cross-region claim must not
## be assembled from runs a human has not signed off. --allow-unlocked exists
## so the stage can be exercised before the PI records the rows, and it is
## stamped onto every output row so such a table can never be mistaken for a
## citable one.
resolve_run <- function(region) {
    if (!allow_unaccepted) {
        return(require_accepted_upstream(MODULE, opts$cohort, region)$run_id)
    }
    acc <- tryCatch(require_accepted_upstream(MODULE, opts$cohort, region)$run_id,
                    error = function(e) NA_character_)
    if (!is.na(acc)) return(acc)
    cand <- list.files(file.path(module_root, "_m", "runs"),
                       pattern = paste0("^scz-", opts$cohort, "-", region, "-"))
    if (length(cand) == 0) stop("No Module 09 run found for ", region)
    sort(cand, decreasing = TRUE)[[1L]]
}
runs <- vapply(regions, resolve_run, character(1))

read_axis <- function(region) {
    f <- file.path(module_root, "_m", "runs", runs[[region]], "results",
                   "architecture-axis-tests.tsv")
    if (!file.exists(f)) stop("Missing axis table for ", region, ": ", f)
    d <- as.data.table(fread(f))
    ## The covariate-adjusted model is the endpoint: a reviewer's first question
    ## about a rank difference is whether length, CpG count, GC and mappability
    ## explain it. The unadjusted rank test is carried alongside, never instead.
    a <- d[model == "logistic_adjusted" & term == predictor]
    if (nrow(a) != 1L) {
        stop("Expected exactly one logistic_adjusted row for '", predictor,
             "' in ", f, ", found ", nrow(a))
    }
    w <- d[model == "wilcoxon_rank_sum"]
    data.table(
        region = region,
        run_id = runs[[region]],
        n_linked = as.integer(a$n_linked),
        n_background = as.integer(a$n_background),
        adj_log_odds_per_sd = as.numeric(a$estimate),
        adj_se = as.numeric(a$std_error),
        adj_or_per_sd = exp(as.numeric(a$estimate)),
        adj_or_lower = exp(as.numeric(a$estimate) - 1.96 * as.numeric(a$std_error)),
        adj_or_upper = exp(as.numeric(a$estimate) + 1.96 * as.numeric(a$std_error)),
        adj_q = as.numeric(a$qvalue),
        adj_significant = isTRUE(as.logical(a$significant)),
        adj_direction = as.character(a$direction),
        unadj_mean_diff = if (nrow(w) == 1L) as.numeric(w$estimate) else NA_real_,
        unadj_q = if (nrow(w) == 1L) as.numeric(w$qvalue) else NA_real_,
        unadj_direction = if (nrow(w) == 1L) as.character(w$direction) else NA_character_,
        covariates = as.character(a$covariates),
        is_primary_region = isTRUE(as.logical(a$is_primary_region)))
}
axis <- rbindlist(lapply(regions, read_axis), use.names = TRUE)

## ------------------------------------------------------------ upstream freshness
##
## A cross-region decision assembled from per-region runs that cite SUPERSEDED
## upstreams is not a current result, however cleanly the regions agree. The
## sealed 2026-09-09 runs cite the pre-rescore Module 02 runs
## (lgv-AA-{region}-20260823) while the accepted ones are
## lgv-AA-{region}-rescore-20260913 -- the rescore landed four days after them --
## so their axis scores are computed on a score this repo no longer accepts.
## Recorded per region, and citability is withdrawn when any region is stale.
## Acceptance alone is not enough: a run can be accepted and still be built on an
## upstream that has since been superseded.
freshness <- rbindlist(lapply(regions, function(re) {
    man <- fread(file.path(module_root, "_m", "runs", runs[[re]], "manifest.tsv"),
                 colClasses = "character")
    cited <- man$value[man$field == "upstream_local_genetic_variance_run_id"]
    cited <- if (length(cited) == 1L) as.character(cited) else NA_character_
    acc <- tryCatch(
        require_accepted_upstream("02_local_genetic_variance", opts$cohort, re)$run_id,
        error = function(e) NA_character_)
    data.table(region = re, run_id = runs[[re]],
               module_02_cited = cited, module_02_accepted = acc,
               upstream_current = identical(cited, acc))
}))
n_stale <- sum(!freshness$upstream_current)
if (n_stale > 0L) {
    message("[15] WARNING: ", n_stale, " of ", nrow(freshness),
            " region runs cite a superseded Module 02 run:")
    for (i in which(!freshness$upstream_current)) {
        message("  ", freshness$region[i], ": cited ", freshness$module_02_cited[i],
                ", accepted ", freshness$module_02_accepted[i])
    }
    message("  The concordance below is PROVISIONAL. Re-run Module 09 against the ",
            "accepted Module 02 runs before either decision is cited.")
}


## ------------------------------------------------------- donor non-independence
##
## Measured, not asserted. This is the number that decides how the result may be
## described, so it is computed from the catalogs the runs actually cite rather
## than quoted from memory.
donor_list <- function(region) {
    man <- fread(file.path(module_root, "_m", "runs", runs[[region]],
                           "manifest.tsv"), colClasses = "character")
    rid <- man$value[man$field == "upstream_vmr_catalog_run_id"]
    if (length(rid) != 1L || is.na(rid)) return(character())
    f <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs", rid,
                   "vmr", "donors_plink.txt")
    if (!file.exists(f)) return(character())
    d <- fread(f, header = FALSE)
    unique(as.character(d[[ncol(d)]]))
}
donors <- lapply(regions, donor_list)
names(donors) <- regions
overlap <- rbindlist(lapply(utils::combn(regions, 2L, simplify = FALSE), function(p) {
    a <- donors[[p[1]]]; b <- donors[[p[2]]]
    u <- length(union(a, b))
    data.table(region_a = p[1], region_b = p[2],
               n_a = length(a), n_b = length(b),
               n_shared = length(intersect(a, b)),
               jaccard = if (u > 0) length(intersect(a, b)) / u else NA_real_)
}))
n_all_three <- length(Reduce(intersect, donors))
n_union <- length(Reduce(union, donors))

## ------------------------------------- the Module 04 connector, precisely bounded
##
## PI 2026-09-18. Two DIFFERENT statements must not be collapsed:
##
##   (a) SCZ-linked VMRs occupy the lower end of the local-genetic-control axis;
##   (b) SCZ-linked VMRs are directly depleted for a given repeat/repressive
##       annotation.
##
## (b) is only sayable where the DIRECT linked-vs-background annotation test was
## run and is claimable. It is -- 05_integration_enrichment.R compares
## mean_linked against mean_background per annotation -- but claimability is
## narrower than significance: caudate's repeat and H3K9me3 rows carry an
## interpretation_constraint (batch/GC entanglement), so a blanket "depleted of
## repressive chromatin" overstates what may be claimed. This table therefore
## reports per annotation IN HOW MANY REGIONS the depletion is CLAIMABLE, and the
## manuscript wording is read off it rather than written from the p-values.
read_annot <- function(region) {
    f <- file.path(module_root, "_m", "runs", runs[[region]], "results",
                   "integration-enrichment.tsv")
    if (!file.exists(f)) return(data.table())
    d <- as.data.table(fread(f))
    keep <- intersect(c("annotation", "direction", "significant", "claimable",
                        "mean_linked", "mean_background", "qvalue",
                        "interpretation_constraint"), names(d))
    d <- d[, keep, with = FALSE]
    d[, region := region][]
}
annot <- rbindlist(lapply(regions, read_annot), use.names = TRUE, fill = TRUE)
annot_concordance <- if (nrow(annot)) {
    annot[, .(
        n_regions_tested = .N,
        n_regions_depleted_direction = sum(direction == "depleted_in_scz_linked",
                                           na.rm = TRUE),
        n_regions_significant = sum(significant %in% TRUE),
        n_regions_claimable = sum(claimable %in% TRUE),
        claimable_regions = paste(sort(region[claimable %in% TRUE]), collapse = ","),
        constrained_regions = paste(sort(region[nzchar(
            interpretation_constraint %||% "")]), collapse = ","),
        ## The only statement the manuscript may make about this annotation.
        permitted_statement = {
            cr <- sort(region[claimable %in% TRUE])
            if (length(cr) == 0L) {
                "no claimable difference in any region"
            } else {
                paste0("SCZ-linked VMRs are DEPLETED for ", .BY$annotation,
                       " in ", paste(cr, collapse = ", "),
                       " (direct linked-vs-background test; not an enrichment)")
            }
        }
    ), by = annotation][order(-n_regions_claimable, annotation)]
} else data.table()
if (nrow(annot_concordance)) {
    annot_concordance[, blanket_repressive_chromatin_claim_permitted := FALSE]
    annot_concordance[, note := paste0(
        "Claimability is narrower than significance. Do NOT write a blanket ",
        "'depleted of repressive chromatin' statement; cite the annotations and ",
        "regions listed here.")]
}

## ----------------------------------------------------------- the H1 evaluation
dir_ok_each <- axis$adj_direction == dir_required
n_dir_ok <- sum(dir_ok_each, na.rm = TRUE)
## Reported, and only required if the PI re-locks uniformity as a requirement.
direction_consistent <- all(dir_ok_each, na.rm = TRUE)

supported <- axis$adj_significant & dir_ok_each
n_supported <- sum(supported, na.rm = TRUE)
## "Support in at least one region OTHER THAN caudate": the batch confound means
## caudate can contribute but never carry.
supported_outside <- sum(supported & axis$region != support_outside, na.rm = TRUE)

## ------------------------------------------ 1. the region-general association
region_general_met <- n_supported >= min_support && supported_outside >= 1L

## ------------------------------------------- 2. genuine region-specific departures
##
## A departure must be EVIDENCE, not a magnitude gap or a nonsignificant region.
## Two admissible rules, both from config:
##
##   sign_reversal_with_support  the region's effect runs the other way AND is
##                               FDR-significant doing so. Unambiguous.
##   magnitude_difference        leave-one-out: this region against the OTHERS
##                               POOLED, using independent SEs. The regions share
##                               donors, which makes the true SE of a difference
##                               SMALLER than the independent calculation, so this
##                               test UNDERSTATES significance -- conservative for
##                               declaring region-specificity, which is the safe
##                               direction. Stated on the output row so the
##                               conservatism is not mistaken for precision.
dep_rules <- as.character(rule$region_specific_departure_rules)
dep_alpha <- as.numeric(rule$region_specific_departure_alpha)
departure_veto <- isTRUE(rule$departure_vetoes_region_general_claim)

departures <- rbindlist(lapply(seq_len(nrow(axis)), function(i) {
    others <- setdiff(seq_len(nrow(axis)), i)
    wo <- 1 / axis$adj_se[others]^2
    pooled_others <- sum(wo * axis$adj_log_odds_per_sd[others]) / sum(wo)
    se_others <- sqrt(1 / sum(wo))
    d <- axis$adj_log_odds_per_sd[i] - pooled_others
    se_d <- sqrt(axis$adj_se[i]^2 + se_others^2)
    z <- d / se_d
    p <- 2 * stats::pnorm(-abs(z))
    reversal <- !dir_ok_each[i] && isTRUE(axis$adj_significant[i])
    magdiff <- is.finite(p) && p < dep_alpha
    hit <- c(sign_reversal_with_support = reversal,
             magnitude_difference = magdiff)[dep_rules]
    data.table(
        region = axis$region[i],
        log_odds_this_region = axis$adj_log_odds_per_sd[i],
        log_odds_other_regions_pooled = pooled_others,
        difference = d, difference_se_independent = se_d,
        difference_z = z, difference_p = p,
        direction_matches_general = dir_ok_each[i],
        significant_in_this_region = isTRUE(axis$adj_significant[i]),
        rule_sign_reversal_with_support = reversal,
        rule_magnitude_difference = magdiff,
        is_region_specific_departure = any(hit, na.rm = TRUE),
        departure_rule_met = paste(names(hit)[which(hit)], collapse = ","),
        difference_test_is_conservative_shared_donors = TRUE,
        alpha = dep_alpha)
}))
departure_regions <- departures[is_region_specific_departure == TRUE, region]
n_departures <- length(departure_regions)

## H1 -- the region-general claim. Uniformity is not part of it, and a departure
## does not veto it unless the PI locks that.
h1_met <- region_general_met && !(departure_veto && n_departures > 0L)

## Reported for completeness and explicitly NOT a gate: the regions share
## donors, so this interval is anticonservative. It exists so a reader who looks
## for a pooled number finds it labelled rather than computing their own.
w <- 1 / axis$adj_se^2
pooled <- sum(w * axis$adj_log_odds_per_sd) / sum(w)
pooled_se <- sqrt(1 / sum(w))
Q <- sum(w * (axis$adj_log_odds_per_sd - pooled)^2)
Q_df <- nrow(axis) - 1L
Q_p <- stats::pchisq(Q, Q_df, lower.tail = FALSE)

## The reading names WHICH of the two questions was answered and how, so a
## region-general association with a departure cannot be read as either a clean
## uniform result or a failure.
end <- if (identical(dir_required, "lower_in_scz_linked")) "lower_end" else "higher_end"
reading <- if (n_supported < min_support) {
    "no_region_general_association_under_supported"
} else if (supported_outside < 1L) {
    "association_supported_only_in_the_batch_confounded_region"
} else if (n_departures == 0L) {
    paste0("region_general_association_scz_linked_vmrs_at_the_", end,
           "_of_the_local_control_axis_no_region_specific_departure")
} else {
    paste0("region_general_association_scz_linked_vmrs_at_the_", end,
           "_of_the_local_control_axis_with_region_specific_departure_in_",
           paste(departure_regions, collapse = "_"))
}

concordance <- data.table(
    cohort = opts$cohort,
    endpoint = "cross_region_axis_concordance",
    predictor = predictor,
    n_regions = nrow(axis),
    regions = paste(axis$region, collapse = ","),
    direction_required = dir_required,
    n_regions_direction_ok = n_dir_ok,
    ## Described, not required: biological uniformity is not a gate
    ## (PI 2026-09-18).
    direction_consistent_in_all_regions = direction_consistent,
    biological_uniformity_required = FALSE,
    ## The two answers, kept apart.
    region_general_association = region_general_met,
    n_region_specific_departures = n_departures,
    region_specific_departure_regions = paste(departure_regions, collapse = ","),
    departure_vetoes_region_general_claim = departure_veto,
    n_regions_with_statistical_support = n_supported,
    min_regions_with_statistical_support = min_support,
    support_required_outside = support_outside,
    n_supporting_regions_outside = supported_outside,
    h1_met = h1_met,
    reading = reading,
    ## --- non-independence, carried on the row ---
    n_donors_in_all_regions = n_all_three,
    n_donors_union = n_union,
    max_pairwise_jaccard = suppressWarnings(max(overlap$jaccard, na.rm = TRUE)),
    regions_are_independent_replicates = FALSE,
    independence_note = paste0(
        "The regions share donors (", n_all_three, " donors in all ",
        nrow(axis), " regions; union ", n_union,
        "). These are tissues within one cohort, NOT independent replicates; ",
        "report concordance, never a pooled p-value."),
    ## --- descriptive pooled estimate, NOT a gate ---
    pooled_log_odds_fixed_effect = pooled,
    pooled_se_assumes_independence = pooled_se,
    pooled_or_per_sd = exp(pooled),
    pooled_is_descriptive_only = TRUE,
    pooled_se_is_anticonservative_due_to_shared_donors = TRUE,
    heterogeneity_Q = Q, heterogeneity_df = Q_df, heterogeneity_p = Q_p,
    fdr_alpha = alpha,
    accepted_upstreams_only = !allow_unaccepted,
    n_region_runs_on_superseded_module_02 = n_stale,
    upstreams_current_in_all_regions = n_stale == 0L,
    citable = !allow_unaccepted && n_stale == 0L)

## ----------------------------- decision 2, and only decision 2
##
## Decision 1 (the caudate magnitude claim) is NOT recomputed here: it is a
## per-region call that 12_apply_gates.R makes from Module 08 tier 3. It is
## carried alongside purely so the two decisions can be read together, and the
## independence of the two is stated as a field rather than left to the reader.
## NOT YET COMPUTED is not the same as FAILED. A run whose gate stage has not run
## leaves no scz-decision.tsv, and treating that absence as a failure would
## report SUPPLEMENT_OR_OMIT for a run that is simply still in flight.
per_region_ok <- vapply(regions, function(re) {
    f <- file.path(module_root, "_m", "runs", runs[[re]], "results",
                   "scz-decision.tsv")
    if (!file.exists(f)) return(NA)
    d <- as.data.table(fread(f))
    v <- if ("scz_application_retention" %in% names(d)) {
        as.character(d$scz_application_retention[[1]])
    } else if ("main_text_retention" %in% names(d)) {
        ## The 2026-09-09 runs predate the split and recorded the fused field.
        as.character(d$main_text_retention[[1]])
    } else NA_character_
    !identical(v, "SUPPLEMENT_OR_OMIT")
}, logical(1))

caudate_claim <- {
    f <- file.path(module_root, "_m", "runs", runs[[support_outside]], "results",
                   "scz-decision.tsv")
    if (support_outside %in% names(runs) && file.exists(f)) {
        d <- as.data.table(fread(f))
        if ("caudate_magnitude_claim" %in% names(d)) {
            as.character(d$caudate_magnitude_claim[[1]])
        } else "NOT_RECORDED_PRE_SPLIT_RUN"
    } else "NOT_RECORDED_PRE_SPLIT_RUN"
}

n_pending_gates <- sum(is.na(per_region_ok))
retention <- if (n_pending_gates > 0L) {
    ## Distinguish "the per-region gates have not run yet" from "a region failed
    ## its gates". Only the second is a reason to demote the application.
    "PENDING_PER_REGION_GATES"
} else if (any(per_region_ok %in% FALSE)) {
    "SUPPLEMENT_OR_OMIT"
} else if (h1_met) "RETAIN_MAIN_TEXT" else "SUPPLEMENT_OR_OMIT"

decisions <- data.table(
    cohort = opts$cohort,
    decision_1_caudate_magnitude_claim = caudate_claim,
    decision_1_gated_on = "08_region_donor_generalization tier 3 reading",
    decision_2_scz_application_retention = retention,
    decision_2_gated_on = paste0("axis direction in all regions + support in >=",
                                 min_support, " including >=1 outside ",
                                 support_outside),
    decision_2_independent_of_decision_1 = TRUE,
    n_regions_with_gate_not_yet_run = n_pending_gates,
    h1_met = h1_met,
    region_general_association = region_general_met,
    n_region_specific_departures = n_departures,
    region_specific_departure_regions = paste(departure_regions, collapse = ","),
    biological_uniformity_required = FALSE,
    h1_statement = paste0(
        "SCZ-linked VMRs show ", dir_required, " local genetic-control scores ",
        "with statistical support in >= ", min_support,
        " regions including >= 1 outside ", support_outside,
        ". Uniformity across regions is NOT required; region-specific departures ",
        "are tested separately and qualify the scope of the claim rather than ",
        "vetoing it."),
    regions_are_independent_replicates = FALSE,
    upstreams_current_in_all_regions = n_stale == 0L,
    provisional_pending_module_09_rerun = n_stale > 0L,
    citable = !allow_unaccepted && n_stale == 0L)

if (nrow(annot_concordance)) {
    write_atomic(annot_concordance, file.path(out_dir,
        paste0("scz-annotation-cross-region-concordance-", opts$cohort, ".tsv")))
}

write_atomic(departures, file.path(out_dir,
    paste0("scz-axis-region-specific-departures-", opts$cohort, ".tsv")))

write_atomic(freshness, file.path(out_dir,
    paste0("scz-region-upstream-freshness-", opts$cohort, ".tsv")))

write_atomic(axis, file.path(out_dir,
    paste0("scz-axis-per-region-", opts$cohort, ".tsv")))
write_atomic(overlap, file.path(out_dir,
    paste0("scz-region-donor-overlap-", opts$cohort, ".tsv")))
write_atomic(concordance, file.path(out_dir,
    paste0("scz-axis-cross-region-concordance-", opts$cohort, ".tsv")))
write_atomic(decisions, file.path(out_dir,
    paste0("scz-application-decisions-", opts$cohort, ".tsv")))

print(axis[, .(region, n_linked, n_background,
               or_per_sd = round(adj_or_per_sd, 3),
               ci = sprintf("%.3f-%.3f", adj_or_lower, adj_or_upper),
               q = signif(adj_q, 2), adj_significant, adj_direction)])
cat("\n[15] donor overlap -- these are NOT independent replicates:\n")
print(overlap[, .(region_a, region_b, n_a, n_b, n_shared,
                  jaccard = round(jaccard, 2))])
cat(sprintf("  %d donors in all %d regions; union %d\n",
            n_all_three, nrow(axis), n_union))
cat(sprintf("\n[15] Q1 region-general: direction %s in %d/%d; support in %d (need >=%d), %d outside %s -> %s\n",
            dir_required, n_dir_ok, nrow(axis), n_supported, min_support,
            supported_outside, support_outside, region_general_met))
cat(sprintf("[15] Q2 region-specific departures: %d of %d regions%s\n",
            n_departures, nrow(axis),
            if (n_departures) paste0(" (", paste(departure_regions, collapse = ", "), ")") else ""))
print(departures[, .(region, log_odds_this_region,
                     others_pooled = round(log_odds_other_regions_pooled, 4),
                     diff = round(difference, 4), p = signif(difference_p, 3),
                     reversal = rule_sign_reversal_with_support,
                     magdiff = rule_magnitude_difference,
                     departure = is_region_specific_departure)])
cat("  Uniformity is NOT required. A departure qualifies the scope of the\n")
cat("  region-general claim; it does not veto it. The difference test uses\n")
cat("  independent SEs and is therefore CONSERVATIVE under shared donors.\n")
cat(sprintf("[15] H1 met: %s -- %s\n", h1_met, reading))
if (nrow(annot_concordance)) {
    cat("\n[15] Module 04 connector -- claimable in how many regions:\n")
    print(annot_concordance[, .(annotation, n_regions_depleted_direction,
                                n_regions_significant, n_regions_claimable,
                                claimable_regions)])
    cat("  A blanket 'depleted of repressive chromatin' claim is NOT permitted;\n")
    cat("  cite the annotation and the regions where it is claimable.\n")
}
cat(sprintf("\n[15] DECISION 1 caudate magnitude claim  : %s\n", caudate_claim))
cat(sprintf("[15] DECISION 2 scz application retention: %s\n", retention))
cat("  The two are independent: decision 2 may hold while the caudate magnitude\n",
    "  attenuates after n-matching. Attenuation is NOT bias in the caudate\n",
    "  estimate -- Module 08 tier 3 shows donor count is a plausible major\n",
    "  contributor to the larger magnitude, nothing more.\n", sep = "")
if (allow_unaccepted) {
    message("[15] --allow-unlocked: outputs are marked citable=FALSE. ",
            "Acceptance is a human act (AGENTS.md 6).")
}
if (n_stale > 0L) {
    cat(sprintf("\n[15] PROVISIONAL: %d of %d region runs are built on a superseded\n",
                n_stale, nrow(freshness)))
    cat("     Module 02 run, so citable=FALSE regardless of acceptance. Re-run\n")
    cat("     Module 09 against the accepted rescore runs to make this citable.\n")
}
cat("\n[15] wrote ", out_dir, "\n", sep = "")
