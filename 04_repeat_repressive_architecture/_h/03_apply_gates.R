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
survives <- function(d) {
    prim <- d[analysis_set == "primary"]
    if (nrow(prim) != 1 || is.na(prim$p) || prim$p > ALPHA) return(FALSE)
    sens <- d[analysis_set != "primary"]
    if (nrow(sens) == 0) return(FALSE)
    all(!is.na(sens$p) & sens$p <= ALPHA & sign(sens$estimate) == sign(prim$estimate))
}

primary_pred <- "local_snp_contribution_score_z"
outcomes <- c("h3k9me3_frac", "quiescent_frac", "line_l1_frac")

claims <- rbindlist(lapply(outcomes, function(o) {
    per_region <- rbindlist(lapply(regions, function(r) {
        d <- res[outcome == o & predictor == primary_pred & region == r]
        prim <- d[analysis_set == "primary"]
        highmap <- d[analysis_set == "high_mappability"]
        data.table(
            region = r,
            estimate = if (nrow(prim) == 1) prim$estimate else NA_real_,
            p = if (nrow(prim) == 1) prim$p else NA_real_,
            highmap_estimate = if (nrow(highmap) == 1) highmap$estimate else NA_real_,
            highmap_p = if (nrow(highmap) == 1) highmap$p else NA_real_,
            survives = survives(d)
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

out_dir <- file.path(repo_root(), MODULE, "_m", "runs", run_ids[1], "results")
write_atomic(claims, file.path(out_dir, "interpretation-claims.tsv"))
write_atomic(res, file.path(out_dir, "association-results-all-regions.tsv"))

print(claims)

## Overlap alone never implies activity, expression, or retrotransposition
## (AGENTS.md 7.4). Emitted with the table so it travels with the numbers.
writeLines(c(
    "Interpretation constraints carried by this table:",
    "  - Outcomes are genomic OVERLAP. Overlap does not imply repeat activity,",
    "    expression, or retrotransposition competence.",
    "  - Estimates are per 1 SD of the within-cell local SNP contribution score",
    "    among eligible, within-domain loci; they are not PVE differences.",
    "  - Claims exceeding `permitted_claim` are not supported by this analysis."
), file.path(out_dir, "interpretation-constraints.txt"))
