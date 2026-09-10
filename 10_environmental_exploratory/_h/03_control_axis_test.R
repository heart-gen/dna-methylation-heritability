#!/usr/bin/env Rscript
#### 10_environmental_exploratory -- exposure association along the control axis ####
##
## Usage:
##   Rscript _h/03_control_axis_test.R --run-id env-AA-caudate-YYYYMMDD
##
## Stage B, and the reason this module can exist at all.
##
## v1 asked this question with 01.fishers_enrichment.py: Fisher's exact tests of
## exposure-associated VMRs against the Heritable / Non-heritable / Low
## prediction categories. Those categories come from `h2_unscaled` and
## `r_squared_cv` (AGENTS.md 3, retired) and encode the genetically-anchored vs
## exposure-associated binary (AGENTS.md 2.3, banned). The question survives the
## grouping variable: it becomes a position on Module 02's continuous
## within-cell rank.
##
## Three models, all prespecified in config/environmental.yml.
##
##   1. PRIMARY, threshold-free: -log10(p) ~ score_z + technical covariates.
##      Continuous on both sides, so the headline statement does not depend on
##      where the FDR cut lands.
##   2. Wilcoxon of score_z, FDR-significant vs not.
##   3. Logistic on the same indicator with the same covariates.
##
## (2) and (3) are secondary because they need a threshold; they are kept
## because a reviewer will ask for the grouped form, and because they are the
## nearest legal analogue of what v1 reported.
##
## Everything here is within ONE cohort x region cell. The score is a within-cell
## midrank percentile, uniform by construction, and AGENTS.md 7.6 forbids
## comparing it across cells at all -- so there is no pooled model and no
## cross-region contrast, by design rather than by omission.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "10_environmental_exploratory", "_h")),
    "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "10_environmental_exploratory"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mf("cohort"); region <- mf("region")
env <- load_run_config("environmental", run_dir)
predictor <- env$testing$architecture_predictor
axis_covs <- as.character(unlist(env$testing$axis_covariates))

assoc <- fread(file.path(run_dir, "results", "vmr-exposure-association.tsv"))

score <- load_local_genetic_control(mf("upstream_local_genetic_variance_run_id"),
                                    region = region, cohort = cohort)
if (!predictor %in% names(score)) {
    stop("Module 02 table has no column '", predictor, "'")
}

## Technical covariates come from Module 04's feature matrix, which is the one
## place they are computed for this VMR set. Joining rather than recomputing
## keeps GC content and mappability identical to the values Module 04 published.
feat_run <- mf("upstream_repeat_architecture_run_id")
feat_f <- file.path(repo_root(), "04_repeat_repressive_architecture", "_m",
                    "runs", feat_run, "results", "vmr-features.tsv")
if (!file.exists(feat_f)) stop("Missing Module 04 features: ", feat_f)
feat <- fread(feat_f)
banned <- intersect(unlist(env$forbidden_columns), c(names(feat), names(score)))
if (length(banned)) {
    stop("Upstream table carries superseded estimator column(s): ",
         paste(banned, collapse = ", "), " (AGENTS.md 3)")
}
for (tbl in list(score, feat)) {
    audit <- intersect(unlist(env$audit_only_columns), names(tbl))
    ## AGENTS.md 7.2 keeps the raw estimate for audit inside Module 02. Dropping
    ## it here means it cannot reach a model or a reported number.
    if (length(audit)) tbl[, (audit) := NULL]
}
missing_covs <- setdiff(axis_covs, names(feat))
if (length(missing_covs)) {
    stop("Module 04 feature table lacks declared axis covariate(s): ",
         paste(missing_covs, collapse = ", "))
}

dt <- merge(assoc, score[, c("vmr_id", predictor), with = FALSE],
            by = "vmr_id", all.x = TRUE)
dt <- merge(dt, feat[, c("vmr_id", axis_covs), with = FALSE],
            by = "vmr_id", all.x = TRUE)

## Denominators are explicit (AGENTS.md 11). A VMR dropped for want of a
## covariate is counted, not quietly absent.
dt[, complete := stats::complete.cases(.SD),
   .SDcols = c(predictor, axis_covs, "p_joint")]

## -log10(p) is heavy-tailed and p can be exactly 0 at machine precision, which
## would make the response infinite. Clamp at the smallest representable double
## and record that it happened.
eps <- .Machine$double.xmin
dt[, y := -log10(pmax(p_joint, eps))]

covar_rhs <- paste(axis_covs, collapse = " + ")

test_one <- function(ex, st) {
    d <- dt[exposure == ex & stratum == st & complete == TRUE]
    n_total <- nrow(dt[exposure == ex & stratum == st])
    n_used <- nrow(d)
    n_sig <- sum(d$significant, na.rm = TRUE)
    base <- data.table(
        exposure = ex, stratum = st, cohort = cohort, region = region,
        run_id = opts$run_id,
        predictor = predictor,
        n_vmrs_in_family = n_total, n_vmrs_modelled = n_used,
        n_significant = n_sig,
        n_dropped_incomplete = n_total - n_used
    )
    if (n_used < 100L) {
        return(cbind(base, data.table(status = "too_few_complete_vmrs")))
    }

    ## 1. PRIMARY -- continuous, threshold-free.
    f1 <- stats::as.formula(paste("y ~", predictor, "+", covar_rhs))
    m1 <- stats::lm(f1, data = d)
    c1 <- summary(m1)$coefficients[predictor, ]

    ## 2. Rank test. Reported whatever n_significant is, but flagged when the
    ##    significant set is too small for the comparison to mean anything.
    w <- if (n_sig >= 10L && n_sig <= n_used - 10L) {
        stats::wilcox.test(d[[predictor]][d$significant],
                           d[[predictor]][!d$significant], exact = FALSE)
    } else NULL

    ## 3. Logistic, same covariates.
    g <- if (n_sig >= 10L && n_sig <= n_used - 10L) {
        f3 <- stats::as.formula(paste("significant ~", predictor, "+", covar_rhs))
        fit <- try(stats::glm(f3, data = d, family = stats::binomial()), silent = TRUE)
        if (inherits(fit, "try-error")) NULL else fit
    } else NULL
    g_co <- if (!is.null(g)) summary(g)$coefficients[predictor, ] else NULL

    cbind(base, data.table(
        status = "ok",
        primary_beta = unname(c1["Estimate"]),
        primary_se = unname(c1["Std. Error"]),
        primary_p = unname(c1["Pr(>|t|)"]),
        primary_model = paste(deparse(f1), collapse = " "),
        wilcoxon_p = if (!is.null(w)) unname(w$p.value) else NA_real_,
        wilcoxon_median_score_significant =
            if (!is.null(w)) stats::median(d[[predictor]][d$significant]) else NA_real_,
        wilcoxon_median_score_other =
            if (!is.null(w)) stats::median(d[[predictor]][!d$significant]) else NA_real_,
        logistic_log_or = if (!is.null(g_co)) unname(g_co["Estimate"]) else NA_real_,
        logistic_se = if (!is.null(g_co)) unname(g_co["Std. Error"]) else NA_real_,
        logistic_p = if (!is.null(g_co)) unname(g_co[4]) else NA_real_,
        grouped_tests_skipped_reason =
            if (is.null(w)) "fewer_than_10_vmrs_in_a_group" else NA_character_
    ))
}

pairs <- unique(assoc[, .(exposure, stratum)])
out <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    test_one(pairs$exposure[i], pairs$stratum[i])
}), fill = TRUE)

## The secondary models get their own BH correction across exposures; the
## primary does too, separately. They are different families and are not pooled.
out[status == "ok", primary_fdr := stats::p.adjust(primary_p, method = "BH")]
out[status == "ok" & !is.na(wilcoxon_p),
    wilcoxon_fdr := stats::p.adjust(wilcoxon_p, method = "BH")]
out[status == "ok" & !is.na(logistic_p),
    logistic_fdr := stats::p.adjust(logistic_p, method = "BH")]

out[, `:=`(
    exploratory_supplement_only = TRUE,
    causal_interpretation_allowed = FALSE,
    absolute_pve_interpretation_allowed = FALSE,
    cross_region_comparison_allowed = FALSE,
    collider_flagged = exposure %in% as.character(unlist(
        env$interpretation$collider_flagged_exposures)),
    burden_indicator = exposure %in% as.character(unlist(
        env$interpretation$burden_indicator_exposures))
)]

write_atomic(out, file.path(run_dir, "results", "control-axis-test.tsv"))
append_manifest(list(dir = run_dir), list(
    n_axis_tests = as.character(nrow(out[status == "ok"])),
    axis_predictor = predictor,
    axis_covariates = paste(axis_covs, collapse = ",")
))
print(out[, .(exposure, stratum, n_vmrs_modelled, n_significant,
              primary_beta, primary_p, wilcoxon_p)])
