#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- integration with the 04 architecture ####
##
## Usage:
##   Rscript _h/05_integration_enrichment.R --run-id scz-AA-caudate-YYYYMMDD
##
## AGENTS.md 7.8 requires an explicit integration analysis: are
## schizophrenia-linked VMRs enriched for LINE/L1, H3K9me3, quiescent
## chromatin, high-mappability repeat intervals, and expression or splicing
## coupling? "If this integration is positive, connect the disease application
## directly to the repeat/repressive architecture. If it is null, present Phase
## 7 as a separate proof of disease relevance." Both outcomes are results; this
## stage is not a search for the first one.
##
## The five annotations are a CLOSED, prespecified list from
## config/schizophrenia.yml:integration.test_against, and they form their own
## FDR family. Features are read from Module 04's accepted matrix rather than
## recomputed, so a discrepancy with Module 04's own figures is impossible.
##
## Module 04's acceptance was QUALIFIED: caudate is perfectly confounded with
## sequencing batch, so its LINE/L1 result was withdrawn and its H3K9me3 result
## is GC-entangled. This stage therefore stamps every row with the constraint
## that applies to it, rather than leaving the caveat to the writing.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09_schizophrenia_risk_application", "_h")),
    "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
})

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
scz <- load_run_config("schizophrenia", run_dir)

linkage <- fread(file.path(run_dir, "results", "vmr-scz-linkage.tsv"))
feat_f <- file.path(repo_root(), "04_repeat_repressive_architecture", "_m", "runs",
                    mf("upstream_repeat_architecture_run_id"),
                    "results", "vmr-features.tsv")
if (!file.exists(feat_f)) stop("Missing Module 04 features: ", feat_f)
feat <- fread(feat_f)
dt <- merge(linkage[, .(vmr_id, scz_linked)], feat, by = "vmr_id", all.x = TRUE)
dt[, scz_linked := as.logical(scz_linked)]

## Coupling is the one integration annotation that lives in Module 07 rather
## than Module 04: a VMR is coupled if any enabled modality reported an
## FDR-significant local association for it in the accepted run.
tsc_dir <- file.path(repo_root(), "07_transcription_splicing_coupling", "_m",
                     "runs", mf("upstream_tx_coupling_run_id"), "results")
summ <- list.files(tsc_dir, pattern = "-vmr-summary\\.tsv$", full.names = TRUE)
if (length(summ) == 0) stop("No Module 07 vmr-summary tables in ", tsc_dir)
coupling <- rbindlist(lapply(summ, fread), use.names = TRUE, fill = TRUE)
coupled <- coupling[, .(expression_splicing_coupling =
                            as.numeric(any(as.logical(any_sig_fdr)))),
                    by = vmr_id]
dt <- merge(dt, coupled, by = "vmr_id", all.x = TRUE)
dt[is.na(expression_splicing_coupling), expression_splicing_coupling := 0]

## High-mappability repeat intervals: the AGENTS.md 7.8 annotation is the
## INTERSECTION of "is a repeat interval" and "is confidently mappable". A
## repeat enrichment that lives only in poorly mappable sequence is the
## artefact this annotation exists to separate out.
map_min <- as.numeric(scz$integration$high_mappability_min)
dt[, high_mappability_repeat_intervals :=
       as.numeric(is.finite(line_l1_frac) & line_l1_frac > 0 &
                  is.finite(mappability) & mappability >= map_min)]

feature_cols <- scz$integration$feature_columns
resolve_col <- function(name) {
    if (identical(name, "expression_splicing_coupling")) return(name)
    if (identical(name, "high_mappability_repeat_intervals")) return(name)
    col <- feature_cols[[name]]
    if (is.null(col)) stop("No feature column configured for '", name, "'")
    col
}

## Constraints carried per annotation, from the Module 04 acceptance record.
constraint_for <- function(name) {
    if (!identical(region, "caudate")) return("")
    if (identical(name, "line_l1") ||
        identical(name, "high_mappability_repeat_intervals")) {
        return(paste0("CAUDATE_BATCH_CONFOUNDED: region is perfectly confounded ",
                      "with sequencing batch (AANRI batch 3); Module 04 withdrew ",
                      "its caudate LINE/L1 claim. Not interpretable as regional ",
                      "biology."))
    }
    if (identical(name, "h3k9me3")) {
        return(paste0("CAUDATE_GC_ENTANGLED: Module 04 accepted H3K9me3 in 2/3 ",
                      "regions only, and the caudate estimate is entangled with ",
                      "GC content."))
    }
    ""
}

res <- rbindlist(lapply(unlist(scz$integration$test_against), function(name) {
    col <- resolve_col(name)
    if (!col %in% names(dt)) {
        return(data.table(annotation = name, feature_column = col,
                          model = "not_run", n_linked = NA_integer_,
                          n_background = NA_integer_,
                          mean_linked = NA_real_, mean_background = NA_real_,
                          estimate = NA_real_, statistic = NA_real_,
                          pvalue = NA_real_,
                          not_run_reason = "feature column absent from Module 04 matrix"))
    }
    v <- as.numeric(dt[[col]])
    ok <- is.finite(v)
    a <- v[ok & dt$scz_linked]
    b <- v[ok & !dt$scz_linked]
    if (length(a) < 3 || length(b) < 3 || length(unique(c(a, b))) < 2) {
        return(data.table(annotation = name, feature_column = col,
                          model = "not_run", n_linked = length(a),
                          n_background = length(b),
                          mean_linked = if (length(a)) mean(a) else NA_real_,
                          mean_background = if (length(b)) mean(b) else NA_real_,
                          estimate = NA_real_, statistic = NA_real_,
                          pvalue = NA_real_,
                          not_run_reason = "too few informative VMRs or no variation"))
    }
    w <- suppressWarnings(stats::wilcox.test(a, b, alternative = "two.sided"))
    data.table(annotation = name, feature_column = col,
               model = "wilcoxon_rank_sum",
               n_linked = length(a), n_background = length(b),
               mean_linked = mean(a), mean_background = mean(b),
               estimate = mean(a) - mean(b),
               statistic = unname(w$statistic), pvalue = w$p.value,
               not_run_reason = "")
}), use.names = TRUE, fill = TRUE)

res[, `:=`(run_id = opts$run_id, cohort = cohort, region = region,
           fdr_family = "scz_integration_annotations_per_region",
           high_mappability_min = map_min)]
## The family is the five prespecified annotations of THIS region, corrected
## once. Rows that could not be run carry no p-value and so enter no family.
res[, qvalue := NA_real_]
testable <- res[, is.finite(pvalue)]
if (any(testable)) {
    res[testable, qvalue := p.adjust(pvalue, method = scz$testing$fdr_method)]
}
res[, significant := is.finite(qvalue) & qvalue < as.numeric(scz$testing$fdr_alpha)]
res[, interpretation_constraint := vapply(annotation, constraint_for, character(1))]
## Direction, carried with the estimate. Several of these annotations can come
## out either way, and a significant NEGATIVE difference means SCZ-linked VMRs
## are DEPLETED for the annotation -- the opposite of the enrichment the
## integration analysis was framed to look for. Recording it here is what stops
## "enriched for LINE/L1" being written about a depletion.
res[, direction := fifelse(!is.finite(estimate), NA_character_,
                    fifelse(estimate > 0, "enriched_in_scz_linked",
                            "depleted_in_scz_linked"))]
res[, claim_wording := fcase(
    !(significant %in% TRUE), "no significant difference",
    direction == "enriched_in_scz_linked",
        paste0("SCZ-linked VMRs are ENRICHED for ", annotation),
    default = paste0("SCZ-linked VMRs are DEPLETED for ", annotation,
                     " (this is not an enrichment)"))]
## A constrained annotation may not carry an integration claim even when it is
## statistically significant.
res[, claimable := significant & !nzchar(interpretation_constraint)]
## Whether the AGENTS.md 7.8 framing -- "connect the disease application to the
## repeat/repressive architecture" -- is supported. Only a claimable ENRICHMENT
## does that; a claimable depletion is a different finding and is reported as
## one.
res[, supports_repressive_architecture_link :=
        claimable & direction == "enriched_in_scz_linked"]

setorder(res, qvalue, pvalue, annotation)
write_atomic(res, file.path(run_dir, "results", "integration-enrichment.tsv"))

print(res[, .(annotation, n_linked, estimate, pvalue, qvalue, claimable)])
message("[09] integration: ", sum(res$claimable, na.rm = TRUE), "/", nrow(res),
        " annotations claimable in ", region)
