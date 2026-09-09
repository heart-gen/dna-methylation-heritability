#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- enrichment along the local-control axis ####
##
## Usage:
##   Rscript _h/04_architecture_axis.R --run-id scz-AA-caudate-YYYYMMDD
##
## AGENTS.md 7.8 retention criterion: "enrichment along the new
## local-genetic-control axis in the primary region". The axis is Module 02's
## relative within-cell score. Two things follow from what that score IS:
##
##   - it is a RANK within one cell, so it is compared between VMR groups of the
##     same cell and never across regions (AGENTS.md 7.6 forbids that outright);
##   - it carries no absolute meaning, so the test is "are SCZ-linked VMRs
##     higher on the axis than comparable VMRs", not "do they exceed a PVE".
##
## Two models, both prespecified. The rank test makes no distributional
## assumption; the logistic model asks whether the difference survives the
## technical covariates a reviewer will raise (length, CpG count, GC,
## mappability). Neither is a screen: the comparison is fixed before the run.

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
predictor <- scz$testing$architecture_predictor

linkage <- fread(file.path(run_dir, "results", "vmr-scz-linkage.tsv"))
if (!predictor %in% names(linkage)) {
    stop("Linkage table has no column '", predictor, "'")
}

## Technical covariates come from Module 04's feature matrix, which is the one
## place they are computed for this VMR set. Joining rather than recomputing
## keeps GC content and mappability identical to the values Module 04 published.
feat_run <- mf("upstream_repeat_architecture_run_id")
feat_f <- file.path(repo_root(), "04_repeat_repressive_architecture", "_m",
                    "runs", feat_run, "results", "vmr-features.tsv")
if (!file.exists(feat_f)) stop("Missing Module 04 features: ", feat_f)
feat <- fread(feat_f)
banned <- intersect(unlist(scz$forbidden_columns), names(feat))
if (length(banned)) {
    stop("Module 04 feature table carries superseded estimator column(s): ",
         paste(banned, collapse = ", "))
}

covs <- intersect(unlist(scz$testing$architecture_covariates), names(feat))
dt <- merge(linkage, feat[, c("vmr_id", covs), with = FALSE],
            by = "vmr_id", all.x = TRUE)
dt[, scz_linked := as.logical(scz_linked)]

score <- dt[[predictor]]
ok <- is.finite(score)
res <- list()

## ------------------------------------------------------------- rank test
w <- suppressWarnings(stats::wilcox.test(score[ok & dt$scz_linked],
                                         score[ok & !dt$scz_linked],
                                         alternative = "two.sided",
                                         conf.int = FALSE))
res[[length(res) + 1]] <- data.table(
    model = "wilcoxon_rank_sum",
    term = "scz_linked",
    estimate = mean(score[ok & dt$scz_linked]) -
               mean(score[ok & !dt$scz_linked]),
    estimate_label = "mean score_z difference (linked - background)",
    std_error = NA_real_,
    statistic = unname(w$statistic),
    pvalue = w$p.value,
    n_linked = sum(ok & dt$scz_linked),
    n_background = sum(ok & !dt$scz_linked),
    covariates = "")

## ------------------------------------------------- covariate-adjusted model
## The score is the PREDICTOR and linkage the outcome, deliberately: linkage is
## a fixed positional fact about published loci, so the direction of the model
## cannot be read as methylation causing anything.
fit_dt <- dt[ok & complete.cases(dt[, c(predictor, covs), with = FALSE])]
if (nrow(fit_dt) > 0 && length(unique(fit_dt$scz_linked)) == 2) {
    form <- stats::as.formula(paste("scz_linked ~", predictor,
                                    if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""))
    fit <- stats::glm(form, data = fit_dt, family = stats::binomial())
    co <- summary(fit)$coefficients
    if (predictor %in% rownames(co)) {
        res[[length(res) + 1]] <- data.table(
            model = "logistic_adjusted",
            term = predictor,
            estimate = co[predictor, "Estimate"],
            estimate_label = "log-odds of SCZ linkage per 1 SD of the score",
            std_error = co[predictor, "Std. Error"],
            statistic = co[predictor, "z value"],
            pvalue = co[predictor, "Pr(>|z|)"],
            n_linked = sum(fit_dt$scz_linked),
            n_background = sum(!fit_dt$scz_linked),
            covariates = paste(covs, collapse = ","))
    }
}

out <- rbindlist(res, use.names = TRUE, fill = TRUE)
out[, `:=`(run_id = opts$run_id, cohort = cohort, region = region,
           predictor = predictor,
           is_primary_region = identical(region, scz$primary_region),
           fdr_family = "scz_architecture_axis_per_region")]
out[, qvalue := p.adjust(pvalue, method = scz$testing$fdr_method)]
out[, significant := is.finite(qvalue) & qvalue < as.numeric(scz$testing$fdr_alpha)]

## Direction is recorded, not assumed. The retention criterion is stated in
## config/analysis_thresholds.yml:phase7_scz as "enrichment toward higher
## predictability OR interpretable contrast", so a significant difference in
## either direction satisfies it -- but a DEPLETION must never be reported
## with the word "enrichment", and the only reliable way to prevent that is to
## carry the direction and the permitted wording with the estimate.
out[, direction := fifelse(!is.finite(estimate), NA_character_,
                    fifelse(estimate > 0, "higher_in_scz_linked",
                            "lower_in_scz_linked"))]
out[, contrast_label := fcase(
    !(significant), "no significant contrast",
    direction == "higher_in_scz_linked",
        "SCZ-linked VMRs sit HIGHER on the local-genetic-control axis (enrichment)",
    default =
        "SCZ-linked VMRs sit LOWER on the local-genetic-control axis (DEPLETION, not enrichment)")]
write_atomic(out, file.path(run_dir, "results", "architecture-axis-tests.tsv"))

## Descriptive distribution, so the effect can be read without refitting.
desc <- dt[ok, .(n = .N,
                 mean_score_z = mean(get(predictor)),
                 median_score_z = stats::median(get(predictor)),
                 q25 = stats::quantile(get(predictor), 0.25),
                 q75 = stats::quantile(get(predictor), 0.75)),
           by = .(scz_linked)]
write_atomic(desc, file.path(run_dir, "results", "architecture-axis-descriptive.tsv"))

print(out[, .(model, estimate, pvalue, qvalue, n_linked, n_background)])
message("[09] architecture axis tested in ", region,
        if (identical(region, scz$primary_region)) " (primary region)" else "")
