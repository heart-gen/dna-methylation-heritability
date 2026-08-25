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

## The adjustment set is DRIVEN BY CONFIG, not restated here. A hard-coded list
## silently diverges from config/repeat_annotations.yml the moment a covariate is
## added, and the model that gets fitted stops being the model that was locked.
## `realized` maps each declared covariate name to the column 01_build_features.R
## actually produced, and to the model term for it.
REALIZED_TERM <- list(
    vmr_length = "log(vmr_length)", cpg_count = "cpg_count",
    cpg_density = "cpg_density", gc_content = "gc_content",
    mean_methylation = "mean_methylation",
    methylation_variance = "methylation_variance",
    wgbs_coverage = "wgbs_coverage", tested_snp_count = "tested_snp_count",
    snp_proximity = "snp_proximal_frac", mappability = "mappability",
    segdup_overlap = "segdup_frac",
    problematic_region_overlap = "problematic_frac",
    broad_genomic_annotation = "factor(broad_genomic_annotation)",
    cell_composition_pcs = "cell_composition_r2"
)
declared <- unlist(annot$covariates)
unknown <- setdiff(declared, names(REALIZED_TERM))
if (length(unknown)) {
    stop("config declares covariate(s) with no model term here: ",
         paste(unknown, collapse = ", "))
}
COVARIATES <- unlist(REALIZED_TERM[declared])

## A term whose column is constant carries no information and makes the design
## rank-deficient. Drop it, and RECORD that it was dropped -- silently fitting a
## thinner model than the locked one is the failure this guards against.
term_col <- unlist(lapply(REALIZED_TERM[declared], function(t)
    gsub("^(log|factor)\\(|\\)$", "", t)))
constant <- vapply(term_col, function(k) {
    if (!k %in% names(feat)) return(TRUE)
    length(unique(feat[[k]][!is.na(feat[[k]])])) < 2L
}, logical(1))
if (any(constant)) {
    write_atomic(data.table(covariate = names(constant)[constant],
                            reason = "constant_or_absent"),
                 file.path(run_dir, "results", "dropped-covariates.tsv"))
    warning("Dropped constant/absent covariate(s): ",
            paste(names(constant)[constant], collapse = ", "), call. = FALSE)
    COVARIATES <- COVARIATES[!constant]
}

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
    exclude_segdups = function(d) d[segdup_frac == 0],
    ## config sensitivities.exclude_snp_proximal_cpgs. At VMR level the
    ## realizable form is to drop loci whose span is SNP-proximal, since the
    ## outcome is a per-VMR overlap rather than a per-CpG measurement.
    exclude_snp_proximal = function(d) d[snp_proximal_frac == 0],
    ## config sensitivities.cell_composition_rna_music: restrict to loci whose
    ## methylation is not dominated by cell composition. The adjustment is
    ## already in every model; this asks whether the result survives removing
    ## the composition-driven loci entirely.
    low_cell_composition = function(d) {
        thr <- stats::quantile(d$cell_composition_r2, 0.75, na.rm = TRUE)
        d[is.finite(cell_composition_r2) & cell_composition_r2 <= thr]
    }
)

fit_one <- function(d, outcome, predictor, set_name) {
    if (nrow(d) < 50) {
        return(data.table(analysis_set = set_name, outcome = outcome,
                          predictor = predictor, n = nrow(d), n_fitted = 0L,
                          estimate = NA_real_, se = NA_real_, z = NA_real_,
                          p = NA_real_, note = "fewer than 50 loci"))
    }
    binary <- endsWith(outcome, "_any")
    ## config primary_model: quasibinomial logit on the overlap FRACTION, with
    ## the binary any-overlap call as the legacy-comparable secondary form.
    fam <- if (binary) stats::binomial() else
        stats::quasibinomial(link = annot$primary_model$link %||% "logit")
    form <- stats::as.formula(paste(outcome, "~", predictor, "+",
                                    paste(COVARIATES, collapse = " + ")))
    fit <- tryCatch(stats::glm(form, data = d, family = fam),
                    error = function(e) NULL)
    if (is.null(fit)) {
        return(data.table(analysis_set = set_name, outcome = outcome,
                          predictor = predictor, n = nrow(d), n_fitted = 0L,
                          estimate = NA_real_, se = NA_real_, z = NA_real_,
                          p = NA_real_, note = "model did not converge"))
    }
    ## Heteroskedasticity-robust SEs. VMRs are not donors, so there is no donor
    ## clustering here, but the proportion outcomes are strongly
    ## heteroskedastic and model-based SEs would be optimistic.
    ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = annot$primary_model$vcov %||% "HC3"))
    if (!predictor %in% rownames(ct)) {
        return(data.table(analysis_set = set_name, outcome = outcome,
                          predictor = predictor, n = nrow(d), n_fitted = 0L,
                          estimate = NA_real_, se = NA_real_, z = NA_real_,
                          p = NA_real_, note = "predictor dropped (collinear)"))
    }
    ## n is the set size; n_fitted is how many rows the model actually used.
    ## They differ for the chromatin outcomes, which are NA on VMRs that did not
    ## lift to hg19. Reporting only n would overstate the evidence.
    data.table(analysis_set = set_name, outcome = outcome, predictor = predictor,
               n = nrow(d), n_fitted = stats::nobs(fit),
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
