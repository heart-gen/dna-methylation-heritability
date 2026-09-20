#!/usr/bin/env Rscript
#### 08 Stage 05 -- the acceptance gate ####
##
## Usage:
##   Rscript _h/05_apply_gates.R --run-id rdg-AA-crossregion-YYYYMMDD
##
## Evaluates the ten criteria in
## config/region_donor_generalization.yml:gate.criteria and writes a terminal
## decision. Note what the criteria do and do not test: they check that every
## result is correctly TIERED and correctly CONSTRAINED, not that any result is
## positive. A null cross-region replication and a null region difference both
## pass this gate, exactly as Module 06 passed with
## sldsc_supports_brain_enrichment = FALSE. A gate that required a finding would
## be a gate on the conclusion.
##
## Acceptance itself remains a human act: this stage produces the evidence, and
## nothing here writes a README row (AGENTS.md 6).

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
    if (length(v) != 1L) NA_character_ else as.character(v[[1L]])
}
if (run_is_sealed(manifest)) {
    stop("Run is sealed and immutable: ", opts$run_id)
}

cfg <- load_config("region_donor_generalization")
out_dir <- file.path(run_dir, "results")
criteria_expected <- as.character(config_get(cfg, "gate.criteria"))
decision_pass <- as.character(config_get(cfg, "gate.decision_pass"))

read_out <- function(f, required = TRUE) {
    p <- file.path(out_dir, f)
    if (!file.exists(p)) {
        if (required) stop("Gate needs ", f, "; run the earlier stages first")
        return(NULL)
    }
    as.data.table(fread(p))
}

tiers        <- read_out("tiers.tsv")
tests        <- read_out("cross-region-tests.tsv")
replication  <- read_out("cross-region-replication.tsv")
rank_agree   <- read_out("cross-region-rank-agreement.tsv")
cr_summary   <- read_out("cross-region-summary.tsv")
identified   <- read_out("identified-difference.tsv")
id_summary   <- read_out("identified-difference-summary.tsv")
descriptive  <- read_out("descriptive-confounded-regions.tsv")
donor_group  <- read_out("donor-group-concordance.tsv")
tier3_on     <- identical(toupper(mval("tier3_enabled")), "TRUE")
ds_reps      <- read_out("caudate-downsampling-replicates.tsv", required = tier3_on)
ds_summary   <- read_out("caudate-downsampling-summary.tsv", required = tier3_on)

regions <- strsplit(mval("regions"), ",", fixed = TRUE)[[1]]
valid_tiers <- c(tiers$tier, "donor_group_concordance")

## Every table this module publishes must carry a tier. Listed explicitly
## rather than globbed, so adding an output without tiering it fails the gate
## instead of slipping past an overly clever glob.
tiered_outputs <- list(
    "cross-region-tests.tsv" = tests,
    "cross-region-replication.tsv" = replication,
    "cross-region-rank-agreement.tsv" = rank_agree,
    "identified-difference.tsv" = identified,
    "descriptive-confounded-regions.tsv" = descriptive,
    "donor-group-concordance.tsv" = donor_group)
if (tier3_on) tiered_outputs[["caudate-downsampling-replicates.tsv"]] <- ds_reps

check <- function(name, value, detail) {
    data.table(criterion = name, pass = isTRUE(value), detail = detail)
}

results <- list()

## 1 -------------------------------------------------------------------------
bad_tier <- vapply(names(tiered_outputs), function(nm) {
    d <- tiered_outputs[[nm]]
    if (is.null(d) || !nrow(d)) return(TRUE)
    if (!"tier" %in% names(d)) return(TRUE)
    u <- unique(as.character(d$tier))
    length(u) != 1L || !u %in% valid_tiers
}, logical(1))
results[[1]] <- check(
    "every_output_carries_exactly_one_tier", !any(bad_tier),
    if (any(bad_tier)) paste("untiered or multi-tiered:",
                             paste(names(tiered_outputs)[bad_tier],
                                   collapse = ", "))
    else paste(length(tiered_outputs), "outputs, each exactly one known tier"))

## 2 -------------------------------------------------------------------------
present <- sort(unique(as.character(tests$region)))
results[[2]] <- check(
    "cross_region_replication_computed_for_all_three_regions",
    setequal(present, regions) && nrow(replication) > 0,
    paste0("regions present: ", paste(present, collapse = ","),
           "; expected: ", paste(regions, collapse = ","),
           "; ", nrow(replication), " tests assessed"))

## 3 -------------------------------------------------------------------------
contrast <- as.character(config_get(cfg, "identified_difference.contrast"))
confounded <- as.character(
    config_get(cfg, "interpretation.technically_confounded_regions.all_outcomes"))
claimed_contrast <- unique(as.character(identified$contrast))
results[[3]] <- check(
    "identified_difference_restricted_to_dlpfc_hippocampus",
    length(claimed_contrast) == 1L &&
        identical(claimed_contrast, paste(contrast[[1]], "minus", contrast[[2]])) &&
        !any(confounded %in% contrast) &&
        identified[difference_claimed == TRUE & !is.finite(delta_z), .N] == 0L,
    paste0("contrast: ", paste(claimed_contrast, collapse = "/"),
           "; confounded regions excluded: ",
           paste(confounded, collapse = ",")))

## 4 -------------------------------------------------------------------------
if (tier3_on) {
    want <- as.integer(mval("tier3_replicates"))
    got <- nrow(ds_reps)
    results[[4]] <- check(
        "caudate_downsampling_replicates_complete",
        got == want && all(is.finite(ds_reps$retained_fraction)) &&
            isTRUE(ds_summary$caudate_remains_batch_confounded[[1]]) &&
            !isTRUE(ds_summary$residual_excess_may_be_called_biological[[1]]),
        paste0(got, " of ", want, " replicates; reading '",
               ds_summary$reading[[1]], "'; batch-confounding flag retained"))
} else {
    results[[4]] <- check("caudate_downsampling_replicates_complete", FALSE,
                          "tier 3 is disabled; the gate cannot pass without it")
}

## 5 -------------------------------------------------------------------------
policy <- donor_group_inference_policy()
results[[5]] <- check(
    "donor_group_policy_is_concordance_only",
    identical(policy$mode, "concordance_only") &&
        !policy$ancestry_effect_claim_allowed &&
        !policy$cross_group_raw_score_comparison &&
        nrow(donor_group) > 0 &&
        all(as.character(donor_group$donor_group_inference) == "concordance_only") &&
        !any(toupper(as.character(donor_group$ancestry_effect_claim_allowed)) == "TRUE"),
    paste0("policy '", policy$mode, "'; ancestry_effect_claim_allowed=",
           policy$ancestry_effect_claim_allowed, "; ", nrow(donor_group),
           " regions carry it"))

## 6 -------------------------------------------------------------------------
flag_true <- function(d, col) {
    !is.null(d) && col %in% names(d) &&
        any(toupper(as.character(d[[col]])) == "TRUE")
}
results[[6]] <- check(
    "no_pooled_rank_or_pooled_r2_emitted",
    !flag_true(donor_group, "pooled_rank_emitted") &&
        !flag_true(donor_group, "pooled_r2_emitted"),
    "no pooled rank and no pooled r2 in any donor-group output")

## 7 -------------------------------------------------------------------------
## Tier 1 must report rank and direction, not levels. The guard is that the
## replication table carries no raw cross-region score column and that the
## summary says so.
banned <- grep("^(score|pve|local_snp_contribution)", names(replication),
               value = TRUE)
results[[7]] <- check(
    "no_cross_region_raw_score_comparison",
    length(banned) == 0L &&
        !flag_true(cr_summary, "raw_score_comparison_emitted") &&
        identical(as.character(cr_summary$agreement_basis[[1]]),
                  "rank_and_direction"),
    paste0("agreement_basis=", cr_summary$agreement_basis[[1]],
           if (length(banned)) paste("; banned columns:",
                                     paste(banned, collapse = ",")) else ""))

## 8 -------------------------------------------------------------------------
## Set aside, not dropped: the confounded region's estimates must still be
## present and must still be flagged as supporting nothing.
results[[8]] <- check(
    "confounded_caudate_reported_separately_not_dropped",
    nrow(descriptive) > 0 &&
        setequal(unique(as.character(descriptive$region)), confounded) &&
        !any(toupper(as.character(descriptive$supports_any_claim)) == "TRUE") &&
        identified[, .N] > 0,
    paste0(nrow(descriptive), " ", paste(confounded, collapse = "/"),
           " rows retained, none supporting a claim"))

## 9 -------------------------------------------------------------------------
results[[9]] <- check(
    "reliability_ceiling_attached_to_every_concordance_figure",
    "reliability_ceiling" %in% names(donor_group) &&
        all(is.finite(as.numeric(donor_group$reliability_ceiling))) &&
        all(is.finite(as.numeric(donor_group$fraction_of_ceiling))) &&
        all(as.numeric(donor_group$spearman_score) <=
                as.numeric(donor_group$reliability_ceiling) + 1e-6),
    paste0("ceilings present for ", nrow(donor_group),
           " regions; observed agreement at ",
           paste(sprintf("%.1f%%", 100 * as.numeric(donor_group$fraction_of_ceiling)),
                 collapse = "/"), " of ceiling"))

## 10 ------------------------------------------------------------------------
## Criterion 2 only asks that every region is PRESENT. It passed a run whose
## completeness flag was FALSE for all 154 tests, so tier 1 reported zero
## replication and the gate still read 9/9 (6c2a930e4). A null replication is
## a finding and must pass; a replication assessed over zero complete tests is
## no finding at all. Completeness is re-derived here from each test's own
## region count rather than trusted from the flag stage 01 wrote, so a defect
## in that flag is caught even when the counts it feeds agree with it.
## Hoisted: `replication` carries a `regions` column, so length(regions)
## inside replication[...] would count ROWS -- the exact defect this
## criterion exists to catch.
n_regions_expected <- length(regions)
n_complete_rederived <- replication[as.integer(n_regions) == n_regions_expected, .N]
n_complete_flagged <- replication[as.logical(complete_across_regions) %in% TRUE, .N]
n_complete_summary <- as.integer(cr_summary$n_tests_complete_all_sets[[1]])
n_claim_complete <- as.integer(cr_summary$n_claim_complete[[1]])
results[[10]] <- check(
    "cross_region_completeness_nonvacuous",
    n_complete_rederived > 0L &&
        n_complete_flagged == n_complete_rederived &&
        identical(n_complete_summary, n_complete_rederived) &&
        isTRUE(n_claim_complete > 0L),
    paste0(n_complete_rederived, " of ", nrow(replication),
           " tests observed in all ", n_regions_expected, " regions (flag: ",
           n_complete_flagged, "; summary: ", n_complete_summary, "); ",
           n_claim_complete, " claim-family tests complete"))

qc <- rbindlist(results, use.names = TRUE)

## The criteria evaluated must be exactly the criteria the config commissioned.
## A criterion silently dropped from this stage would make the gate weaker than
## the config claims, which is the failure mode a gate cannot have.
if (!setequal(qc$criterion, criteria_expected)) {
    stop("Gate criteria do not match config/region_donor_generalization.yml.\n",
         "  missing here: ",
         paste(setdiff(criteria_expected, qc$criterion), collapse = ", "),
         "\n  extra here:   ",
         paste(setdiff(qc$criterion, criteria_expected), collapse = ", "))
}
qc <- qc[match(criteria_expected, criterion)]

decision <- if (all(qc$pass)) decision_pass else "FAIL_REGION_DONOR_GENERALIZATION_QC"
qc[, decision := decision]
write_atomic(qc, file.path(out_dir, "region-donor-generalization-qc.tsv"))

decision_dt <- data.table(
    run_id = opts$run_id,
    cohort = mval("cohort"),
    regions = mval("regions"),
    n_criteria = nrow(qc),
    n_passed = sum(qc$pass),
    decision = decision,
    ## Said out loud: this gate certifies TIERING and CONSTRAINT, not a result.
    gate_certifies = "tiering and interpretation constraints, not a positive finding",
    accepted_by = NA_character_,
    accepted_on = NA_character_)
write_atomic(decision_dt,
             file.path(out_dir, "region-donor-generalization-decision.tsv"))

print(qc[, .(criterion, pass, detail)])
cat("\n")
print(decision_dt[, .(n_passed, n_criteria, decision)])
if (!all(qc$pass)) {
    message("[gate] FAILED: ",
            paste(qc[pass == FALSE, criterion], collapse = ", "))
} else {
    message("[gate] ", decision, " -- ", nrow(qc), "/", nrow(qc),
            " criteria pass.")
}
message("  Acceptance is a human act. Nothing is citable until the PI records ",
        "this run in 08_region_donor_generalization/README.md (AGENTS.md 6).")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
