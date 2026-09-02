#!/usr/bin/env Rscript
#### 04_repeat_repressive_architecture -- cross-region interpretation gates ####
##
## Usage:
##   Rscript _h/03_apply_gates.R --cohort AA --run-ids rra-AA-caudate-...,rra-AA-dlpfc-...,rra-AA-hippocampus-...
##
## This script converts model output into the sentences the manuscript is
## allowed to write, and it is the reason 04 has a step after the models.
## AGENTS.md 7.4 sets the rules; encoding them here means the claim is derived,
## not asserted:
##
##   - H3K9me3 and quiescent enrichment may be called "shared across regions"
##     only if DIRECTION AND INFERENCE survive the locked sensitivities in ALL
##     THREE regions.
##   - LINE/L1 may be called multi-region if it survives in AT LEAST TWO. If it
##     survives only in caudate, the permitted claim is "caudate-specific".
##   - A DLPFC reversal or null after the high-mappability restriction is never
##     hidden; it is written into the output table as its own field.
##   - The concentration claim needs BOTH arms: enrichment in a repressive
##     outcome AND depletion in accessible chromatin. Without the second arm the
##     permitted noun is "enriched in <track>", not "concentrated in repressive
##     chromatin", because enrichment alone is equally consistent with VMRs
##     sitting in gene deserts where every repressive track is wide.
##   - "Constitutive heterochromatin" requires the primary repressive signal to
##     exceed the Polycomb/bivalent control. If they are comparable, the
##     permitted noun is "repressed chromatin".
##
## The output is a claims table with an explicit `permitted_claim` string per
## outcome. If a result does not clear its gate, the permitted claim says so.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "04_repeat_repressive_architecture"

opts <- parse_v2_args(require = c("cohort", "run_ids"))
annot <- load_config("repeat_annotations")
## From config, not restated: the interpretation gates in
## config/repeat_annotations.yml are the locked decision rule.
ALPHA <- as.numeric(annot$interpretation$alpha %||% 0.05)
## Family members are gated on the BH-corrected q; controls, which are outside
## the family and therefore uncorrected, are gated on raw p. Previously this
## script tested raw p for everything, silently discarding the correction
## 02_test_association.R had just computed.
GATE_FAMILY  <- annot$multiple_testing$gate_statistic$family %||% "q"
GATE_CONTROL <- annot$multiple_testing$gate_statistic$outside_family %||% "p"

run_ids <- trimws(strsplit(opts$run_ids, ",", fixed = TRUE)[[1]])
res <- rbindlist(lapply(run_ids, function(rid) {
    f <- file.path(repo_root(), MODULE, "_m", "runs", rid, "results",
                   "association-results.tsv")
    if (!file.exists(f)) stop("No association results for run ", rid)
    d <- fread(f); d[, run_id := rid]; d
}))

regions <- sort(unique(res$region))
if (length(regions) != 3) {
    stop("Interpretation gates are defined over all three regions; got: ",
         paste(regions, collapse = ", "),
         "\n  Run all three cells before applying the gates.")
}

#' Does this outcome survive in this region?
#'
#' "Survives" means the primary estimate is significant AND every locked
#' sensitivity keeps the same sign and remains significant. A sensitivity that
#' flips the sign is a failure even if it is significant -- that is a
#' contradiction, not a replication.
#' `stat` is the column the PRIMARY estimate is judged on -- q for family
#' members, p for controls. Sensitivities are always judged on p: they are
#' re-fits of a hypothesis already counted in the family, so correcting them
#' again would penalize the same test twice.
survives <- function(d, stat = GATE_FAMILY) {
    prim <- d[analysis_set == "primary"]
    if (nrow(prim) != 1) return(FALSE)
    pv <- prim[[stat]]
    if (is.na(pv) || pv > ALPHA) return(FALSE)
    sens <- d[analysis_set != "primary"]
    if (nrow(sens) == 0) return(FALSE)
    all(!is.na(sens$p) & sens$p <= ALPHA & sign(sens$estimate) == sign(prim$estimate))
}

#' One-sided p for a prespecified direction. A control with the wrong sign gets
#' no credit, which is the point of declaring the direction in advance.
one_sided_p <- function(estimate, p, direction) {
    if (is.na(estimate) || is.na(p)) return(NA_real_)
    right_sign <- if (identical(direction, "negative")) estimate < 0 else estimate > 0
    if (right_sign) p / 2 else 1 - p / 2
}

primary_pred <- "local_snp_contribution_score_z"
## The gated outcomes are exactly the BH family. Controls are evaluated
## separately below, under their own rules.
outcomes <- unlist(annot$multiple_testing$family)
CONTROLS <- names(annot$multiple_testing$outside_family)

claims <- rbindlist(lapply(outcomes, function(o) {
    per_region <- rbindlist(lapply(regions, function(r) {
        d <- res[outcome == o & predictor == primary_pred & region == r]
        prim <- d[analysis_set == "primary"]
        highmap <- d[analysis_set == "high_mappability"]
        data.table(
            region = r,
            estimate = if (nrow(prim) == 1) prim$estimate else NA_real_,
            p = if (nrow(prim) == 1) prim$p else NA_real_,
            q = if (nrow(prim) == 1) prim$q else NA_real_,
            highmap_estimate = if (nrow(highmap) == 1) highmap$estimate else NA_real_,
            highmap_p = if (nrow(highmap) == 1) highmap$p else NA_real_,
            survives = survives(d, GATE_FAMILY)
        )
    }))

    n_surv <- sum(per_region$survives)
    surv_regions <- per_region$region[per_region$survives]

    required <- if (o == "line_l1_frac") {
        annot$interpretation$line_l1_multiregion_requires
    } else {
        annot$interpretation$shared_requires_regions
    }
    required <- as.integer(required %||% 3L)

    claim <- if (n_surv >= required && n_surv == length(regions)) {
        paste0("shared across all three regions (n=", n_surv, ")")
    } else if (n_surv >= required) {
        paste0("multi-region: ", paste(surv_regions, collapse = ", "))
    } else if (n_surv == 1 && identical(surv_regions, "caudate")) {
        "caudate-specific"
    } else if (n_surv == 0) {
        "not supported: no region survives the locked sensitivities"
    } else {
        paste0("below the gate (", n_surv, "/", required,
               " regions); may be described only as suggestive in ",
               paste(surv_regions, collapse = ", "))
    }

    ## The DLPFC-specific disclosure, surfaced as its own field so it cannot be
    ## dropped by summarizing the table.
    dl <- per_region[region == "dlpfc"]
    dlpfc_note <- if (nrow(dl) == 1 && !is.na(dl$estimate) && !is.na(dl$highmap_estimate)) {
        if (sign(dl$highmap_estimate) != sign(dl$estimate)) {
            "DLPFC direction REVERSES after the high-mappability restriction; must be reported"
        } else if (!is.na(dl$highmap_p) && dl$highmap_p > ALPHA) {
            "DLPFC becomes null after the high-mappability restriction; must be reported"
        } else NA_character_
    } else NA_character_

    cbind(data.table(outcome = o, regions_surviving = n_surv,
                     regions_required = required,
                     permitted_claim = claim, dlpfc_disclosure = dlpfc_note),
          data.table(t(stats::setNames(per_region$survives, per_region$region))))
}))

## ------------------------------------------------------ prespecified controls
##
## Outside the BH family (config multiple_testing.outside_family), so gated on
## raw p and, where a direction was declared in advance, on the ONE-SIDED p.
## They do not carry claims of their own; they qualify the claims above.
controls <- rbindlist(lapply(CONTROLS, function(o) {
    spec <- annot$multiple_testing$outside_family[[o]]
    dir_expected <- spec$expected_direction
    rbindlist(lapply(regions, function(r) {
        d <- res[outcome == o & predictor == primary_pred & region == r]
        prim <- d[analysis_set == "primary"]
        if (nrow(prim) != 1) {
            return(data.table(outcome = o, role = spec$role, region = r,
                              estimate = NA_real_, p = NA_real_,
                              p_directional = NA_real_, expected_direction =
                                  dir_expected %||% NA_character_,
                              supports_expectation = FALSE))
        }
        pd <- if (is.null(dir_expected)) prim$p else
            one_sided_p(prim$estimate, prim$p, dir_expected)
        data.table(
            outcome = o, role = spec$role, region = r,
            estimate = prim$estimate, p = prim$p, p_directional = pd,
            expected_direction = dir_expected %||% NA_character_,
            supports_expectation = !is.na(pd) && pd <= ALPHA &&
                (is.null(dir_expected) || survives(d, GATE_CONTROL))
        )
    }))
}))

## Arm two of the concentration claim. Depletion in accessible chromatin must
## hold, in the prespecified direction, in the required number of regions.
acc_cfg <- annot$interpretation$accessible_contrast
acc_required <- as.integer(acc_cfg$regions_required %||% length(regions))
acc <- controls[role == "complementary_contrast" & outcome == "accessible_frac"]
n_acc <- sum(acc$supports_expectation)
concentration_supported <- !isTRUE(acc_cfg$required_for_concentration_claim) ||
    n_acc >= acc_required

## Constitutive vs. Polycomb. Compared on the standardized estimate rather than
## on p, because a larger sample on one track would otherwise decide it.
##
## SIGN-AWARE (config amendment 2026-09-02). Specificity is a contrast, not two
## marginal fits. What defeats "constitutive heterochromatin" is Polycomb
## behaving like the primary: SAME sign and comparable magnitude, which is
## repression generally. An OPPOSITE-signed Polycomb estimate is evidence FOR
## constitutive specificity -- H3K9me3 and H3K27me3 occupy anti-correlated
## compartments -- and must never block the claim. The previous abs() comparison
## was sign-blind and inverted exactly that case.
spec_cfg <- annot$interpretation$specificity_control
POLYCOMB <- unlist(spec_cfg$polycomb_outcomes %||% "h3k27me3_frac")

signed_mean <- function(dt) {
    v <- dt[!is.na(estimate), estimate]
    if (length(v) == 0) NA_real_ else mean(v)
}
polycomb_signed <- signed_mean(controls[outcome %in% POLYCOMB])
repressive_signed <- signed_mean(
    res[outcome %in% c("h3k9me3_frac", "quiescent_frac") &
        predictor == primary_pred & analysis_set == "primary"])

## Same sign is a precondition for the control to bite at all.
same_sign <- is.finite(polycomb_signed) && is.finite(repressive_signed) &&
    sign(polycomb_signed) == sign(repressive_signed)

constitutive_supported <- !isTRUE(
    spec_cfg$constitutive_claim_requires_stronger_than_polycomb) ||
    !is.finite(polycomb_signed) || !is.finite(repressive_signed) ||
    !same_sign ||
    abs(repressive_signed) > abs(polycomb_signed)

## Reported so the reason is auditable from the log, not inferred from the noun.
constitutive_basis <- if (!is.finite(polycomb_signed)) {
    "no Polycomb estimate available"
} else if (!same_sign) {
    sprintf(paste("Polycomb %.3f is OPPOSITE in sign to repressive %.3f:",
                  "consistent with constitutive specificity"),
            polycomb_signed, repressive_signed)
} else {
    sprintf(paste("Polycomb %.3f and repressive %.3f share sign;",
                  "|repressive| %s |Polycomb|"),
            polycomb_signed, repressive_signed,
            if (abs(repressive_signed) > abs(polycomb_signed)) ">" else "<=")
}

## Fold both qualifiers into the permitted claim, so a reader of the claims
## table alone cannot pick up the strong noun without the evidence for it.
claims[, concentration_supported := concentration_supported]
claims[, constitutive_supported := constitutive_supported]
claims[regions_surviving > 0 & !concentration_supported,
       permitted_claim := paste0(permitted_claim,
           "; ENRICHMENT ONLY -- accessible-chromatin depletion not established (",
           n_acc, "/", acc_required,
           " regions), so \"concentrated in repressive chromatin\" is not supported")]
claims[outcome %in% c("h3k9me3_frac", "quiescent_frac") &
       regions_surviving > 0 & !constitutive_supported,
       permitted_claim := paste0(permitted_claim,
           "; say \"repressed chromatin\", NOT \"constitutive heterochromatin\"",
           " -- ", constitutive_basis)]

out_dir <- file.path(repo_root(), MODULE, "_m", "runs", run_ids[1], "results")
write_atomic(claims, file.path(out_dir, "interpretation-claims.tsv"))
write_atomic(res, file.path(out_dir, "association-results-all-regions.tsv"))
write_atomic(controls, file.path(out_dir, "control-outcomes.tsv"))

print(claims)
print(controls)

## Overlap alone never implies activity, expression, or retrotransposition
## (AGENTS.md 7.4). Emitted with the table so it travels with the numbers.
writeLines(c(
    "Interpretation constraints carried by this table:",
    "  - Outcomes are genomic OVERLAP. Overlap does not imply repeat activity,",
    "    expression, or retrotransposition competence.",
    "  - Estimates are per 1 SD of the within-cell local SNP contribution score",
    "    among eligible, within-domain loci; they are not PVE differences.",
    "  - Claims exceeding `permitted_claim` are not supported by this analysis.",
    "  - Family outcomes are gated on the BH-corrected q; the prespecified",
    "    controls in control-outcomes.tsv are outside that family and are",
    "    reported with raw (directional where declared) p.",
    sprintf("  - Accessible-chromatin depletion: %d/%d regions. Concentration claim %s.",
            n_acc, acc_required,
            if (concentration_supported) "permitted" else "NOT permitted"),
    sprintf("  - Repressive |estimate| %.3f vs Polycomb/bivalent %.3f. %s.",
            repressive_strength, polycomb_strength,
            if (constitutive_supported) "\"Constitutive heterochromatin\" permitted"
            else "Say \"repressed chromatin\"")
), file.path(out_dir, "interpretation-constraints.txt"))
