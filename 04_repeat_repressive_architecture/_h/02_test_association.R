#!/usr/bin/env Rscript
#### 04_repeat_repressive_architecture -- primary and sensitivity models ####
##
## Usage:
##   Rscript _h/02_test_association.R --run-id rra-AA-caudate-20260817
##
## One model family, fitted to the primary analysis set and then re-fitted under
## every locked sensitivity in config/repeat_annotations.yml. The sensitivities
## are not optional add-ons: AGENTS.md 7.4 states enrichment that does not
## survive them is not reported as a result, so they are fitted in the same pass
## and land in the same table. There is no codepath that produces the primary
## estimate alone.
##
## Outcome scale: the fraction outcomes are proportions bounded in [0, 1] with
## mass at both ends, so quasibinomial with a logit link is the default and the
## binary any-overlap outcomes use logistic regression. Both give
## overdispersion-tolerant inference; neither assumes a normal error.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(sandwich)
    library(lmtest)
})

MODULE <- "04_repeat_repressive_architecture"

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
feat <- fread(file.path(run_dir, "results", "vmr-features.tsv"))
annot <- load_config("repeat_annotations")

COVARIATES <- c("log(vmr_length)", "cpg_count", "cpg_density",
                "mean_methylation", "methylation_variance",
                "tested_snp_count", "mappability", "segdup_frac")

OUTCOMES <- c("h3k9me3_frac", "quiescent_frac", "line_l1_frac",
              "h3k9me3_any", "quiescent_any", "line_l1_any")

PREDICTORS <- c(primary = "local_snp_contribution_score_z",
                secondary = if ("r2_pred_oof_z" %in% names(feat)) "r2_pred_oof_z")
PREDICTORS <- PREDICTORS[!vapply(PREDICTORS, is.null, logical(1))]

#' The locked analysis sets. Each returns a subset of `feat`.
analysis_sets <- list(
    primary = function(d) d,
    high_mappability = function(d) {
        d[mappability >= annot$sensitivities$high_mappability$min_mappability]
    },
    exclude_segdups = function(d) d[segdup_frac == 0]
)

fit_one <- function(d, outcome, predictor, set_name) {
    if (nrow(d) < 50) {
        return(data.table(analysis_set = set_name, outcome = outcome,
                          predictor = predictor, n = nrow(d),
                          estimate = NA_real_, se = NA_real_, z = NA_real_,
                          p = NA_real_, note = "fewer than 50 loci"))
    }
    binary <- endsWith(outcome, "_any")
    fam <- if (binary) stats::binomial() else stats::quasibinomial()
    form <- stats::as.formula(paste(outcome, "~", predictor, "+",
                                    paste(COVARIATES, collapse = " + ")))
    fit <- tryCatch(stats::glm(form, data = d, family = fam),
                    error = function(e) NULL)
    if (is.null(fit)) {
        return(data.table(analysis_set = set_name, outcome = outcome,
                          predictor = predictor, n = nrow(d),
                          estimate = NA_real_, se = NA_real_, z = NA_real_,
                          p = NA_real_, note = "model did not converge"))
    }
    ## Heteroskedasticity-robust SEs. VMRs are not donors, so there is no donor
    ## clustering here, but the proportion outcomes are strongly
    ## heteroskedastic and model-based SEs would be optimistic.
    ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
    if (!predictor %in% rownames(ct)) {
        return(data.table(analysis_set = set_name, outcome = outcome,
                          predictor = predictor, n = nrow(d),
                          estimate = NA_real_, se = NA_real_, z = NA_real_,
                          p = NA_real_, note = "predictor dropped (collinear)"))
    }
    data.table(analysis_set = set_name, outcome = outcome, predictor = predictor,
               n = nrow(d),
               estimate = ct[predictor, 1], se = ct[predictor, 2],
               z = ct[predictor, 3], p = ct[predictor, 4], note = NA_character_)
}

results <- rbindlist(lapply(names(analysis_sets), function(set_name) {
    d <- analysis_sets[[set_name]](copy(feat))
    rbindlist(lapply(OUTCOMES, function(o) {
        rbindlist(lapply(PREDICTORS, function(p) fit_one(d, o, p, set_name)))
    }))
}))

## Multiple testing is applied within the primary analysis set across the three
## outcome families, not across the sensitivities -- a sensitivity is a re-fit of
## the same hypothesis, not a new one.
results[, q := NA_real_]
results[analysis_set == "primary" & predictor == "local_snp_contribution_score_z",
        q := stats::p.adjust(p, method = "BH")]

results[, `:=`(region = feat$region[1], population = feat$population[1],
               vmr_set_id = feat$vmr_set_id[1])]

write_atomic(results, file.path(run_dir, "results", "association-results.tsv"))
print(results[analysis_set == "primary"])
message("[04] fitted ", nrow(results), " models across ",
        length(analysis_sets), " analysis sets")
