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
## Arms may add terms the primary does not carry (adjust_cell_composition),
## so the guard runs over the UNION of base and arm terms -- otherwise an
## arm-only term could be constant and enter the design unchecked.
ARM_TERMS <- c(cell_composition_r2 = "cell_composition_r2")
term_col <- c(unlist(lapply(REALIZED_TERM[declared], function(t)
    gsub("^(log|factor)\\(|\\)$", "", t))), ARM_TERMS)
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
    dropped <- names(constant)[constant]
    COVARIATES <- COVARIATES[!names(COVARIATES) %in% dropped]
    ARM_TERMS <- ARM_TERMS[!names(ARM_TERMS) %in% dropped]
}

## The BH family is declared in config, not here, so it cannot drift between
## runs. Everything else fitted below is either a second scale of a family
## member (`_any`) or a prespecified control, and is reported with raw p.
BH_FAMILY <- unlist(annot$multiple_testing$family)
if (length(BH_FAMILY) == 0) stop("multiple_testing.family is empty in config")
CONTROLS <- names(annot$multiple_testing$outside_family)

OUTCOMES <- c(BH_FAMILY, CONTROLS,
              sub("_frac$", "_any", c(BH_FAMILY, CONTROLS)))
missing <- setdiff(OUTCOMES, names(feat))
if (length(missing) > 0) {
    stop("Declared outcomes absent from the feature table: ",
         paste(missing, collapse = ", "),
         "\n  Re-run 01_build_features.R after a config change.")
}

PREDICTORS <- c(primary = "local_snp_contribution_score_z",
                secondary = if ("r2_pred_oof_z" %in% names(feat)) "r2_pred_oof_z")
PREDICTORS <- PREDICTORS[!vapply(PREDICTORS, is.null, logical(1))]

#' The locked analysis sets. Each is a row subset PLUS, optionally, extra
#' adjustment terms: an arm can differ from the primary in what it fits as well
#' as in which loci it fits on.
analysis_sets <- list(
    primary = list(subset = function(d) d, extra_terms = character(0)),
    high_mappability = list(
        subset = function(d) {
            d[mappability >= annot$sensitivities$high_mappability$min_mappability]
        },
        extra_terms = character(0)),
    exclude_segdups = list(subset = function(d) d[segdup_frac == 0],
                           extra_terms = character(0)),
    ## config sensitivities.cell_composition_rna_music: restrict to loci whose
    ## methylation is not dominated by cell composition. This is the SUBSET
    ## form of the composition check; adjust_cell_composition below is the
    ## adjustment form. They are independent and both are gating.
    low_cell_composition = list(
        subset = function(d) {
            thr <- stats::quantile(d$cell_composition_r2, 0.75, na.rm = TRUE)
            d[is.finite(cell_composition_r2) & cell_composition_r2 <= thr]
        },
        extra_terms = character(0)),
    ## config sensitivities.adjust_cell_composition. cell_composition_pcs left
    ## the primary adjustment set in the 2026-09-02 amendment because its
    ## confounder-vs-mediator status is arguable; rather than settle that
    ## silently inside the primary, this arm refits every model with it added.
    adjust_cell_composition = list(subset = function(d) d,
                                   extra_terms = unname(ARM_TERMS))
)

## config sensitivities.exclude_snp_proximal_cpgs, RETIRED from the gating
## conjunction 2026-09-02 as redundant: the identical BED already enters every
## model as the `snp_proximity` covariate (snp_proximal_frac). It is also not a
## faithful realization of the key, which names a CpG-level exclusion --
## snp_proximal_frac is a base-pair overlap fraction of the VMR SPAN, so `== 0`
## selects short VMRs mechanically while every outcome is a length-dependent
## overlap fraction. Still fitted, written to its own file so it can never
## re-enter survives() in 03_apply_gates.R.
DESCRIPTIVE_SETS <- list(
    exclude_snp_proximal = list(subset = function(d) d[snp_proximal_frac == 0],
                                extra_terms = character(0))
)

fit_one <- function(d, outcome, predictor, set_name,
                    terms = COVARIATES) {
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
                                    paste(terms, collapse = " + ")))
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

fit_sets <- function(sets) {
    rbindlist(lapply(names(sets), function(set_name) {
        arm <- sets[[set_name]]
        d <- arm$subset(copy(feat))
        terms <- c(COVARIATES, arm$extra_terms)
        rbindlist(lapply(OUTCOMES, function(o) {
            rbindlist(lapply(PREDICTORS, function(p)
                fit_one(d, o, p, set_name, terms)))
        }))
    }))
}

results <- fit_sets(analysis_sets)

## Multiple testing is applied within the primary analysis set, over the
## declared family only. Three exclusions, each deliberate:
##   - sensitivities, because a sensitivity is a re-fit of the same hypothesis;
##   - the `_any` scale, because it is the same hypothesis measured differently
##     (correcting over both scales would make the family six for three
##     questions, and the two scales are strongly correlated, so the penalty is
##     conservative without being principled);
##   - the prespecified controls, which are directional checks on the
##     interpretation rather than members of the discovery family. Keeping them
##     out is what lets a control be added later without revising any q already
##     reported.
results[, q := NA_real_]
results[analysis_set == "primary" &
        predictor == "local_snp_contribution_score_z" &
        outcome %in% BH_FAMILY,
        q := stats::p.adjust(p, method = "BH")]

## Carried into the results table so a reader of the TSV alone can tell which
## rows were corrected and what each control was expected to do.
results[, `:=`(outcome_role = "secondary_scale", expected_direction = NA_character_)]
results[outcome %in% BH_FAMILY, outcome_role := "bh_family"]
for (o in CONTROLS) {
    spec <- annot$multiple_testing$outside_family[[o]]
    results[outcome == o, outcome_role := spec$role %||% "control"]
    if (!is.null(spec$expected_direction)) {
        results[outcome == o, expected_direction := spec$expected_direction]
    }
}

results[, `:=`(region = feat$region[1], population = feat$population[1],
               vmr_set_id = feat$vmr_set_id[1])]

write_atomic(results, file.path(run_dir, "results", "association-results.tsv"))

## The retired arm, refit and reported but deliberately NOT in
## association-results.tsv: 03_apply_gates.R reads that file and conjoins every
## non-primary analysis set, so a separate file is what keeps this non-gating.
descriptive <- fit_sets(DESCRIPTIVE_SETS)
descriptive[, `:=`(region = feat$region[1], population = feat$population[1],
                   vmr_set_id = feat$vmr_set_id[1], gating = FALSE)]
write_atomic(descriptive,
             file.path(run_dir, "results", "descriptive-snp-proximal.tsv"))
print(results[analysis_set == "primary"])
message("[04] fitted ", nrow(results), " models across ",
        length(analysis_sets), " analysis sets")
