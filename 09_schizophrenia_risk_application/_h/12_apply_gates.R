#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- acceptance gate and retention decision ####
##
## Usage:
##   Rscript _h/12_apply_gates.R --run-id scz-AA-caudate-YYYYMMDD
##
## Two separate questions, kept separate on purpose.
##
## 1. THE GATE: was the analysis conducted correctly? Loci defined without
##    reference to methylation, every stage ran, the tested universe is large
##    enough to interpret, the FDR family stayed its own, no banned quantity
##    entered a model, and colocalization ran on enough ancestry-matched regions
##    with enough shared variants to mean something. A NULL result is not a gate
##    failure -- AGENTS.md 7.8 explicitly provides for a null Phase 7 ("present
##    Phase 7 as a separate proof of disease relevance"), and Modules 06 and 07
##    set the precedent that a reportable null seals.
##
## 2. MAIN-TEXT RETENTION: should the application appear in the main text? That
##    is the five prespecified criteria in config, and it is NOT the gate. One
##    of them -- caudate_not_sample_size_artifact -- depends on the downsampling
##    arm that Module 08 owns. While gates.require_module_08_downsampling is
##    false it is recorded as PENDING_MODULE_08: neither a pass nor a failure, an
##    explicitly unevaluable criterion, and main_text_retention stays
##    PENDING_MODULE_08 no matter how the other four resolve, so the open
##    dependency cannot be lost in the writing. Once the PI flips that flag the
##    criterion is read off the accepted Module 08 run's tier-3 reading, pinned in
##    this run's manifest -- never off the flag itself.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09_schizophrenia_risk_application", "_h")),
    "run_config.R"))

suppressPackageStartupMessages(library(data.table))

MODULE <- "09_schizophrenia_risk_application"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mf("cohort"); region <- mf("region")
smoke <- identical(mf("smoke_run"), "TRUE")
scz <- load_run_config("schizophrenia", run_dir)
gates <- scz$gates
crit <- scz$retention_criteria
is_primary <- identical(region, scz$primary_region)

read_or_empty <- function(name) {
    f <- file.path(run_dir, "results", name)
    if (file.exists(f)) fread(f) else data.table()
}

universe <- read_or_empty("tested-universe.tsv")
rv_summary <- read_or_empty("risk-variant-test-summary.tsv")
axis <- read_or_empty("architecture-axis-tests.tsv")
integration <- read_or_empty("integration-enrichment.tsv")
coupling <- read_or_empty("locus-coupling.tsv")
gtex <- read_or_empty("gtex-support.tsv")
coloc <- fread(file.path(run_dir, "results", "coloc-abf.tsv.gz"))
prio <- read_or_empty("prioritized-loci.tsv")
locus_def <- read_or_empty("locus-definition-summary.tsv")

fail <- character()
num <- function(dt, col, default = 0) {
    if (nrow(dt) == 0 || !col %in% names(dt)) return(default)
    v <- suppressWarnings(as.numeric(dt[[col]][1]))
    if (is.na(v)) default else v
}

## ------------------------------------------------------------- conduct checks
## The design invariant first: if locus definition ever consulted methylation,
## nothing downstream means anything.
if (nrow(locus_def) == 0 || !isTRUE(as.logical(locus_def$methylation_used[1]) == FALSE)) {
    fail <- c(fail, "LOCUS_DEFINITION_NOT_METHYLATION_INDEPENDENT")
}
for (stage_file in c("scz-loci.tsv", "locus-vmr-links.tsv",
                     "risk-variant-cpg-tests.tsv.gz",
                     "architecture-axis-tests.tsv", "integration-enrichment.tsv",
                     "locus-coupling.tsv", "gtex-support.tsv",
                     "coloc-abf.tsv.gz", "prioritized-loci.tsv")) {
    if (!file.exists(file.path(run_dir, "results", stage_file))) {
        fail <- c(fail, paste0("STAGE_OUTPUT_MISSING(", stage_file, ")"))
    }
}
if (!identical(rv_summary$fdr_family[1] %||% "", scz$testing$fdr_family)) {
    fail <- c(fail, "RISK_VARIANT_FDR_FAMILY_NOT_AS_CONFIGURED")
}

if (!smoke) {
    if (num(locus_def, "n_loci_published") < gates$min_loci_tested) {
        fail <- c(fail, "TOO_FEW_LOCI_TESTED")
    }
    if (num(universe, "n_vmrs_linked") < gates$min_vmrs_linked) {
        fail <- c(fail, "TOO_FEW_VMRS_LINKED")
    }
    if (num(rv_summary, "n_pairs_tested") < gates$min_pairs_tested) {
        fail <- c(fail, "TOO_FEW_PAIRS_TESTED")
    }
}

## Every reported test must have produced a finite estimate. A table of NAs
## would otherwise pass silently as "no enrichment".
if (nrow(axis) == 0 || all(!is.finite(axis$estimate))) {
    fail <- c(fail, "NO_ARCHITECTURE_TEST_PRODUCED_AN_ESTIMATE")
}

## ---------------------------------------------------------------- coloc gate
## Evaluated on the ancestry-matched arms ONLY. The cross-ancestry arm can
## neither pass nor fail this gate; including it would let a violated LD
## assumption satisfy the condition that exists to enforce that assumption.
coloc_gate_arms <- if (nrow(coloc)) {
    coloc[gate_eligible == TRUE & ld_ancestry_matched == TRUE]
} else data.table()
n_coloc_evaluated <- if (nrow(coloc_gate_arms))
    coloc_gate_arms[status == "OK", uniqueN(locus_id)] else 0L
coloc_gate_passed <- TRUE
if (!smoke) {
    if (n_coloc_evaluated < gates$min_loci_coloc_tested) {
        fail <- c(fail, "TOO_FEW_LOCI_COLOCALIZATION_EVALUATED")
        coloc_gate_passed <- FALSE
    }
    if (nrow(coloc_gate_arms) == 0) {
        fail <- c(fail, "NO_ANCESTRY_MATCHED_COLOCALIZATION_ARM_RAN")
        coloc_gate_passed <- FALSE
    }
}
## A run whose claimable colocalizations came from an arm that is not
## ancestry-matched would be a contradiction in the output; check it rather
## than trusting stage 09 to have got it right.
if (nrow(coloc) &&
    any(coloc$claimable_colocalization %in% TRUE &
        !(coloc$ld_ancestry_matched %in% TRUE))) {
    fail <- c(fail, "CLAIMABLE_COLOC_FROM_ANCESTRY_UNMATCHED_ARM")
    coloc_gate_passed <- FALSE
}

## ------------------------------------------------------- retention criteria
n_loci_meqtl <- num(rv_summary, "n_loci_with_significant_meqtl")
c1 <- if (n_loci_meqtl >= crit$min_loci_with_cpg_meqtl_support) "PASS" else "FAIL"

axis_sig <- nrow(axis) > 0 && any(axis$significant %in% TRUE)
## config/analysis_thresholds.yml:phase7_scz states this criterion as
## "enrichment_toward_higher_predictability_or_interpretable_contrast", so a
## significant contrast in EITHER direction satisfies it. The direction is
## reported in the status itself, because a depletion that satisfies the
## criterion must not be written up as an enrichment.
axis_dir <- if (axis_sig) {
    sig <- axis[significant %in% TRUE]
    sig$direction[which.max(abs(sig$estimate))]
} else NA_character_
## Evaluated in EVERY region, not only the primary one. PI 2026-09-18: the axis
## contrast is the shared, cross-region finding this module's framing question
## rests on, and scoring it NOT_APPLICABLE outside caudate meant the strongest
## consistent result was assessed only in the one region that is perfectly
## confounded with sequencing batch and so can never carry a claim alone
## (AGENTS.md 8.1). `is_primary` is still recorded, but it no longer decides
## whether the criterion is read.
c2 <- if (!axis_sig) "FAIL" else
      if (identical(axis_dir, "higher_in_scz_linked")) "PASS_ENRICHMENT" else
      "PASS_INTERPRETABLE_CONTRAST_DEPLETION"

n_loci_coupled <- if (nrow(coupling) && "has_transcriptional_coupling" %in% names(coupling))
    sum(coupling$has_transcriptional_coupling %in% TRUE) else 0L
c3 <- if (n_loci_coupled >= crit$min_loci_with_transcriptional_coupling) "PASS" else "FAIL"

## Module 08 owns the caudate downsampling arm. This criterion asks whether
## caudate's excess survives matching its donor count to DLPFC's -- so it is read
## OFF Module 08's tier-3 reading, never inferred from the config flag. The flag
## only says whether the dependency is required yet; it is not evidence, and a
## flag that passed the criterion by itself would let the whole question be
## settled by editing a boolean.
##
## Tier 3's two permitted readings map directly, and the mapping is INVERTED
## relative to the criterion's name: if donor count is a plausible major
## contributor, the caudate excess IS largely a sample-size effect, so
## `caudate_not_sample_size_artifact` FAILS. A third reading is forbidden by
## Module 08's own constraints -- a surviving excess never becomes biological,
## because caudate stays confounded with sequencing batch (AGENTS.md 8.1).
##
## SCOPE (PI 2026-09-18): caudate only. The question is about CAUDATE's
## magnitude, so it is not evaluable for a DLPFC or hippocampus run and must not
## gate one. Before this it was recorded as PENDING_MODULE_08 in all three
## regions, so a caudate sub-question held the entire application open in every
## region.
c4_scope <- as.character(config_get(
    scz, "independent_decisions.caudate_magnitude_claim.scope_regions"))
c4 <- if (!region %in% c4_scope) {
    "NOT_APPLICABLE_NON_CAUDATE_REGION"
} else if (!isTRUE(gates$require_module_08_downsampling)) {
    "PENDING_MODULE_08"
} else {
    m08 <- mf("region_donor_generalization_run_id")
    if (is.na(m08) || !nzchar(m08)) {
        stop("gates.require_module_08_downsampling is true but this run pins no ",
             "Module 08 run. Re-open the run so 00_new_run.R records ",
             "region_donor_generalization_run_id from the accepted 08 run; ",
             "criterion 4 is not evaluable from the flag alone.")
    }
    f08 <- file.path(repo_root(), "08_region_donor_generalization", "_m", "runs",
                     m08, "results", "caudate-downsampling-summary.tsv")
    if (!file.exists(f08)) {
        stop("Module 08 run ", m08, " has no caudate-downsampling-summary.tsv: ",
             f08, "\n  Criterion 4 cannot be evaluated without the tier-3 result.")
    }
    t3 <- fread(f08)
    if (!"reading" %in% names(t3) || nrow(t3) != 1L) {
        stop("Module 08 tier-3 summary must be one row with a `reading` column: ",
             f08)
    }
    switch(as.character(t3$reading[[1]]),
           donor_count_does_not_explain_the_excess = "PASS",
           donor_count_is_a_plausible_major_contributor = "FAIL_SAMPLE_SIZE_ARTIFACT",
           indeterminate_replicates_disagree = "INDETERMINATE_MODULE_08",
           stop("Unrecognized Module 08 tier-3 reading '", t3$reading[[1]],
                "' in ", f08, ". Refusing to map an unknown reading onto this ",
                "criterion."))
}

n_loci_gtex <- if (nrow(gtex) && "has_external_genetic_support" %in% names(gtex)) {
    prio_ids <- if (nrow(prio)) prio[prioritized == TRUE, locus_id] else character()
    sum(gtex$has_external_genetic_support %in% TRUE &
        gtex$locus_id %in% prio_ids)
} else 0L
c5 <- if (n_loci_gtex >= crit$min_loci_with_external_genetic_support) "PASS" else "FAIL"

criteria <- c(loci_with_cpg_meqtl_support = c1,
              enrichment_along_local_control_axis = c2,
              loci_with_transcriptional_coupling = c3,
              caudate_not_sample_size_artifact = c4,
              loci_with_external_genetic_support = c5)

PASSING <- c("PASS", "PASS_ENRICHMENT", "PASS_INTERPRETABLE_CONTRAST_DEPLETION",
             "NOT_APPLICABLE_NON_PRIMARY_REGION",
             "NOT_APPLICABLE_NON_CAUDATE_REGION")

## ------------------------------------------------- TWO INDEPENDENT DECISIONS
##
## PI 2026-09-18. The module's framing question -- does regional variation in
## local genetic control of methylation intersect schizophrenia-relevant
## regulatory biology, and is that relationship shared or region-dependent --
## has two separable answers, and Module 08 tier 3 speaks to only one of them.
## Fusing them into a single main_text_retention let a caudate sub-question
## decide the fate of the whole application in all three regions.
##
## DECISION 1 -- the caudate-specific magnitude claim. Caudate only, gated on
## Module 08 tier 3 and on nothing else.
caudate_magnitude_claim <- if (!region %in% c4_scope) {
    "NOT_APPLICABLE_NON_CAUDATE_REGION"
} else switch(c4,
    PASS = "CAUDATE_MAGNITUDE_CLAIM_SUPPORTED",
    FAIL_SAMPLE_SIZE_ARTIFACT = "CAUDATE_MAGNITUDE_CLAIM_NOT_SUPPORTED",
    PENDING_MODULE_08 = "PENDING_MODULE_08",
    INDETERMINATE_MODULE_08 = "PENDING_MODULE_08_INDETERMINATE",
    stop("Unhandled criterion-4 state for the caudate magnitude claim: ", c4))

## DECISION 2 -- the schizophrenia application. Deliberately EXCLUDES criterion
## 4: the application may stand on a disorder-related pattern shared across the
## three regions even if the caudate magnitude attenuates after n-matching.
## What it cannot be decided on is this run alone -- H1 is a statement about all
## three regions (direction in 3/3, statistical support in >=2 including >=1
## non-caudate), so the cross-region part is resolved by
## _h/15_cross_region_axis_concordance.R over the three sealed runs and recorded
## here as PENDING_CROSS_REGION. A single region can only fail this decision on
## its own per-region criteria, never pass it.
application_criteria <- criteria[names(criteria) != "caudate_not_sample_size_artifact"]
scz_application_retention <- if (!all(application_criteria %in% PASSING)) {
    "SUPPLEMENT_OR_OMIT"
} else "PENDING_CROSS_REGION"

## Retained for continuity with the sealed 2026-09-09 runs and anything that
## reads this field: it now tracks DECISION 2 only, and the caudate magnitude
## claim is reported in its own column rather than folded in here.
main_text_retention <- scz_application_retention

## Module 06's accepted S-LDSC result is null, and
## config/analysis_thresholds.yml:phase7_scz lists
## adds_nothing_beyond_nonsignificant_sldsc as an omit_or_supplement_if
## condition. Carry the upstream finding rather than restating it from memory.
sldsc_f <- file.path(repo_root(), "06_partitioned_heritability", "_m", "runs",
                     mf("upstream_partitioned_h2_run_id"), "results",
                     "partitioned-h2-decision.tsv")
sldsc_supports <- NA
if (file.exists(sldsc_f)) {
    sl <- fread(sldsc_f)
    if ("sldsc_supports_brain_enrichment" %in% names(sl)) {
        sldsc_supports <- as.logical(sl$sldsc_supports_brain_enrichment[1])
    }
}

decision <- if (length(fail) > 0) {
    paste0("FAIL_SCZ_APPLICATION_QC:", paste(fail, collapse = ";"))
} else if (smoke) {
    "PASS_SMOKE_ONLY_NOT_ACCEPTABLE"
} else {
    "PASS_SCZ_APPLICATION_QC"
}

dec <- data.table(
    run_id = opts$run_id, cohort = cohort, region = region,
    decision = decision, smoke_run = smoke,
    is_primary_region = is_primary,
    n_loci_published = num(locus_def, "n_loci_published"),
    n_vmrs_linked = num(universe, "n_vmrs_linked"),
    n_pairs_tested = num(rv_summary, "n_pairs_tested"),
    n_pairs_significant = num(rv_summary, "n_pairs_significant"),
    n_loci_with_cpg_meqtl_support = n_loci_meqtl,
    architecture_axis_significant = axis_sig,
    architecture_axis_direction = axis_dir,
    n_integration_annotations_claimable =
        if (nrow(integration)) sum(integration$claimable %in% TRUE) else 0L,
    ## Only a claimable ENRICHMENT licenses the AGENTS.md 7.8 framing that
    ## connects the disease application to the repeat/repressive architecture.
    n_integration_annotations_enriched =
        if (nrow(integration))
            sum(integration$supports_repressive_architecture_link %in% TRUE) else 0L,
    integration_supports_repressive_architecture_link =
        if (nrow(integration))
            any(integration$supports_repressive_architecture_link %in% TRUE) else FALSE,
    n_loci_with_transcriptional_coupling = n_loci_coupled,
    n_prioritized_loci_with_gtex_support = n_loci_gtex,
    n_prioritized_loci = if (nrow(prio)) sum(prio$prioritized %in% TRUE) else 0L,
    n_loci_coloc_evaluated_ancestry_matched = n_coloc_evaluated,
    n_claimable_colocalizations =
        if (nrow(coloc)) sum(coloc$claimable_colocalization %in% TRUE) else 0L,
    n_coloc_regions_cross_ancestry_exploratory =
        if (nrow(coloc)) sum(!(coloc$ld_ancestry_matched %in% TRUE)) else 0L,
    coloc_gate_passed = coloc_gate_passed,
    ## The only flag that licenses the word "colocalization" in the manuscript.
    coloc_claim_permitted = coloc_gate_passed && length(fail) == 0 && !smoke,
    criterion_loci_with_cpg_meqtl_support = c1,
    criterion_enrichment_along_local_control_axis = c2,
    criterion_loci_with_transcriptional_coupling = c3,
    criterion_caudate_not_sample_size_artifact = c4,
    criterion_loci_with_external_genetic_support = c5,
    ## The two independent decisions. main_text_retention is kept as an alias of
    ## scz_application_retention for continuity; read the two named fields.
    caudate_magnitude_claim = caudate_magnitude_claim,
    scz_application_retention = scz_application_retention,
    main_text_retention = main_text_retention,
    module_08_downsampling_available =
        identical(mf("module_08_downsampling_available"), "TRUE"),
    upstream_sldsc_supports_brain_enrichment = sldsc_supports,
    failures = paste(fail, collapse = ";")
)
write_atomic(dec, file.path(run_dir, "results", "scz-decision.tsv"))
write_atomic(data.table(criterion = names(criteria), status = unname(criteria)),
             file.path(run_dir, "results", "retention-criteria.tsv"))

writeLines(c(
    "Interpretation constraints carried by this run:",
    "",
    "  Design invariants",
    "  - PGC3 loci were defined with no reference to any methylation result;",
    "    locus-definition-summary.tsv records methylation_used = FALSE.",
    paste0("  - Risk-variant x CpG tests form their own FDR family (",
           scz$testing$fdr_family, "); they are never pooled with any other."),
    "  - Association is tested against the relative local SNP contribution",
    "    score. Absolute PVE and the legacy predictability metric are banned.",
    "",
    "  Colocalization",
    paste0("  - Claimable ONLY from the ancestry-matched arms (PGC3 European x ",
           "GTEx v11 European brain), and only while coloc_claim_permitted is ",
           "TRUE. This run: ", dec$coloc_claim_permitted[1], "."),
    "  - The CpG meQTL arm is European GWAS against African-American meQTL.",
    "    coloc assumes both studies share an LD structure and they do not, so",
    "    every row of that arm is stamped CROSS_ANCESTRY_LD_UNMATCHED and is",
    "    EXPLORATORY regardless of its PP4. It becomes gate-eligible only when",
    "    an ancestry-matched meQTL panel exists (config: arms.meqtl).",
    "  - A shared causal variant is not mediation and implies no causal",
    "    ordering between methylation and expression.",
    "  - coloc.susie is a sensitivity check on the single-causal-variant",
    "    assumption. It never sets a claim and never overturns coloc.abf.",
    "",
    "  Direction of effect",
    paste0("  - Architecture axis: ",
           if (is.na(axis_dir)) "no significant contrast" else axis_dir,
           if (identical(axis_dir, "lower_in_scz_linked"))
               paste0(". This is a DEPLETION: SCZ-linked VMRs sit LOWER on ",
                      "the local-genetic-control axis than the matched ",
                      "background. It satisfies the prespecified criterion ",
                      "as an interpretable contrast, and must NOT be ",
                      "written up as enrichment.")
           else "."),
    paste0("  - Integration annotations reaching an ENRICHMENT rather than ",
           "merely a significant difference: ",
           if (nrow(integration))
               sum(integration$supports_repressive_architecture_link %in% TRUE)
           else 0L,
           ". Only these support connecting the disease application to the ",
           "repeat/repressive architecture; significant depletions are a ",
           "different finding and carry their own wording in ",
           "integration-enrichment.tsv:claim_wording."),
    "",
    "  Region-specific",
    if (identical(region, "caudate"))
        paste0("  - Caudate is perfectly confounded with sequencing batch ",
               "(AANRI batch 3). Module 04 withdrew its caudate LINE/L1 claim ",
               "and its caudate H3K9me3 estimate is GC-entangled; the ",
               "corresponding integration rows are marked not claimable.")
    else "  - No region-specific confound recorded for this cell.",
    "",
    "  Open dependencies",
    paste0("  - caudate_not_sample_size_artifact is ", c4, ". It scopes to ",
           "caudate ONLY and gates the caudate magnitude claim alone; it does ",
           "NOT gate the schizophrenia application (PI 2026-09-18)."),
    paste0("  - caudate_magnitude_claim = ", caudate_magnitude_claim,
           ". Module 08 tier 3 answers this and nothing else. NOT_SUPPORTED ",
           "means donor count is a plausible major contributor to caudate's ",
           "larger magnitude -- it does NOT mean the caudate estimate is ",
           "biased, and must never be written that way."),
    paste0("  - scz_application_retention = ", scz_application_retention,
           ". PENDING_CROSS_REGION means this run's own criteria are met and ",
           "the decision now rests on cross-region concordance, which a single ",
           "region cannot settle: see _m/combined/ from ",
           "_h/15_cross_region_axis_concordance.R."),
    paste0("  - Upstream S-LDSC (Module 06) supports brain enrichment: ",
           sldsc_supports, ". config/analysis_thresholds.yml lists ",
           "adds_nothing_beyond_nonsignificant_sldsc as an ",
           "omit_or_supplement_if condition."),
    "",
    "  Forbidden claims (config/schizophrenia.yml:claims.forbidden)",
    paste0("  - ", unlist(scz$claims$forbidden))
), file.path(run_dir, "results", "interpretation-constraints.txt"))

print(dec[, .(decision, n_pairs_significant, n_loci_with_cpg_meqtl_support,
              n_claimable_colocalizations, caudate_magnitude_claim,
              scz_application_retention)])
message("[09] decision ", decision,
        "; caudate magnitude claim ", caudate_magnitude_claim,
        "; application retention ", scz_application_retention)
if (startsWith(decision, "FAIL")) quit(status = 1)
