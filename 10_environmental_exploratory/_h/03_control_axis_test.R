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
## Models, all prespecified in config/environmental.yml.
##
##   1. PRIMARY, threshold-free: omega ~ score_z + technical covariates, where
##      omega is the DEBIASED exposure-explained sum of squares per donor,
##      (SS_exposure - df * sigma2_hat) / n, from stage 2. Continuous on both
##      sides, so the headline does not depend on where the FDR cut lands, and
##      unbiased for the exposure's true contribution whatever the SE.
##   2. DESCRIPTIVE: the same model on -log10(p). Retained because it is what
##      this stage reported before 2026-09-19 and a reader will want the
##      comparison, and flagged mechanically_biased_toward_hypothesis = TRUE.
##   3. NON-GATING ARM: the primary outcome and model with
##      `methylation_variance` dropped and nothing else changed. Holding total
##      variance fixed is what activates the variance-budget arithmetic (see
##      config/environmental.yml:interpretation.variance_budget_limitation), and
##      09b excluded total variance from its primary for the same reason. Fitted
##      on the SAME VMR rows and the SAME bootstrap draws as the primary, so the
##      two are paired and their difference is attributable to that one
##      covariate. It cannot promote or demote the primary.
##   4. Wilcoxon of score_z, FDR-significant vs not.
##   5. Logistic on the same indicator with the same covariates.
##
## (3) and (4) are secondary because they need a threshold; they are kept
## because a reviewer will ask for the grouped form, and because they are the
## nearest legal analogue of what v1 reported. In every run so far they have been
## skipped for want of ten VMRs in a group.
##
## WHY THE OUTCOME CHANGED (PI approval 2026-09-19). The per-VMR exposure model
## has no SNP term, so a VMR with strong local SNP control carries its local
## genetic variance in that model's residual. Its SE is inflated, and any outcome
## whose expectation depends on the SE -- |beta|, its rank, |t|, -log10 p --
## couples to the score mechanically. 09b_aging_application hit this first and
## rejected every SE-dependent outcome (AGENTS.md 7.9). Measured here on the
## 2026-09-19 runs: SE^2 rises with score_z at fixed total methylation variance in
## every pooled family (p 1e-13 to 1e-28).
##
## WHY THE INFERENCE CHANGED. This stage used to report the OLS p-value of a
## regression over 9,000-11,000 VMRs, which treats the VMRs as independent of
## each other and of the donor sample. They are neither. The variance is now
## donor-bootstrap variance + delete-one-chromosome block jackknife variance
## (00_shared/axis_inference.R); on the 2026-09-19 runs the jackknife half alone
## was a median 1.28x the OLS SE.
##
## Everything here is within ONE cohort x region cell. The score is a within-cell
## midrank percentile, uniform by construction, and AGENTS.md 7.6 forbids
## comparing it across cells at all -- so there is no pooled model and no
## cross-region contrast, by design rather than by omission.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
H_DIR <- Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "10_environmental_exploratory", "_h"))
source(file.path(H_DIR, "run_config.R"))
source(file.path(H_DIR, "exposure_model.R"))

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
   .SDcols = c(predictor, axis_covs, "p_joint", "omega")]

## ------------------------------------------------------- bootstrap checkpoint
##
## Stage 2 wrote one donor x VMR phenotype matrix per chromosome, with the donor
## covariate and exposure table. The donor ORDER must be identical across
## chromosomes or the cbind would silently mix donors.
ck_dir <- file.path(run_dir, "results", "checkpoint")
ck_files <- sort(list.files(ck_dir, pattern = "^pheno-chr.*\\.rds$",
                            full.names = TRUE))
if (!length(ck_files)) {
    stop("No stage 2 checkpoint under ", ck_dir, "; the donor bootstrap cannot ",
         "run. Re-run stage 2 with the current code.")
}
cks <- lapply(ck_files, readRDS)
donors <- cks[[1]]$donors
for (k in seq_along(cks)) {
    if (!identical(cks[[k]]$donors, donors)) {
        stop("Checkpoint ", basename(ck_files[k]), " carries a different donor ",
             "order from ", basename(ck_files[1]))
    }
}
Y_all <- do.call(cbind, lapply(cks, `[[`, "Y"))
md_all <- cks[[1]]$md
pc_names <- cks[[1]]$pc_names
stopifnot(identical(as.character(md_all$FID), donors))
arm_drops <- as.character(unlist(env$testing$arm_no_methylation_variance_drops))
arm_covs <- setdiff(axis_covs, arm_drops)
if (!identical(env$testing$primary_axis_estimand,
               "proportional_gradient_ratio_functional")) {
    stop("testing.primary_axis_estimand must be proportional_gradient_ratio_functional")
}
if (!isTRUE(env$testing$absolute_axis_sensitivity)) {
    stop("testing.absolute_axis_sensitivity must be true")
}
arm_declared <- "no_methylation_variance" %in%
    as.character(unlist(env$testing$non_gating_axis_arms))
if (arm_declared) {
    if (!length(arm_drops)) {
        stop("testing.non_gating_axis_arms names no_methylation_variance but ",
             "testing.arm_no_methylation_variance_drops is empty")
    }
    missing_drop <- setdiff(arm_drops, axis_covs)
    if (length(missing_drop)) {
        stop("The arm drops covariate(s) that are not in axis_covariates: ",
             paste(missing_drop, collapse = ", "))
    }
}
B <- as.integer(env$testing$n_bootstrap)
if (!is.finite(B) || B < 2L) stop("config/environmental.yml:testing.n_bootstrap missing")
alpha_ci <- as.numeric(env$testing$fdr_alpha)
message("[10] bootstrap checkpoint: ", length(donors), " donors x ",
        ncol(Y_all), " VMRs from ", length(ck_files), " chromosomes; B = ", B)

## -log10(p) is heavy-tailed and p can be exactly 0 at machine precision, which
## would make the response infinite. Clamp at the smallest representable double
## and record that it happened.
eps <- .Machine$double.xmin
dt[, y := -log10(pmax(p_joint, eps))]

covar_rhs <- paste(axis_covs, collapse = " + ")

## One family: the primary (debiased) and descriptive (-log10 p) axis models with
## combined inference, plus the two threshold-based secondaries.
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

    ## ---------------------------------------------- the donor-level refit design
    ## Built by exposure_model.R, the same code stage 2 fitted with, so the
    ## bootstrap resamples the model whose estimate it is putting a SE on.
    keep_don <- exposure_keep_rows(md_all, ex, st, env, phenotype = NULL)
    d_don <- md_all[keep_don, , drop = FALSE]
    des <- exposure_design(d_don, ex, st, env, pc_names)
    X <- des$X
    n_don <- nrow(X)

    ## Axis rows, restricted to VMRs whose phenotype is complete over this
    ## family's donors. A VMR dropped here is counted, never silently short.
    have <- d$vmr_id %in% colnames(Y_all)
    Yf <- Y_all[keep_don, d$vmr_id[have], drop = FALSE]
    finite_vmr <- colSums(!is.finite(Yf)) == 0L
    d <- d[have][finite_vmr]
    Yf <- Yf[, finite_vmr, drop = FALSE]
    n_axis <- nrow(d)
    if (n_axis < 100L) {
        return(cbind(base, data.table(status = "too_few_vmrs_in_checkpoint")))
    }

    ## Self-check: the matrix path must reproduce stage 2's per-VMR omega, which
    ## came from lm() + anova(). If it does not, the bootstrap would be
    ## resampling a different estimator from the one reported.
    fm0 <- fit_term_matrix(X, Yf, des$cols)
    om0 <- debiased_partial_ss(fm0$ss_term, fm0$df_term, fm0$sigma2, n_don)
    dev <- max(abs(om0 - d$omega))
    rel <- dev / max(abs(d$omega))
    if (!is.finite(rel) || rel > 1e-6) {
        stop("Matrix refit does not reproduce stage 2's omega for ", ex, "@", st,
             " (max|diff| = ", signif(dev, 3), ", relative ", signif(rel, 3), ")")
    }
    stopifnot(fm0$df_term == unique(d$df_num))

    ## PI 2026-09-20. The primary estimand is the PROPORTIONAL gradient,
    ## R = beta / mean(omega), estimated as a ratio functional: each bootstrap
    ## draw and each deleted chromosome recomputes the ratio on its own rows.
    ## That is not a cosmetic rescaling. The overall LEVEL of exposure-explained
    ## variance is a nuisance that varies about 5x between chromosomes and 6x
    ## between donor draws, so the absolute gradient does not replicate across
    ## loci -- measured on the sealed -b runs, its block-jackknife SE is 5-10x
    ## its bootstrap SE and all three of this module's findings go to p ~ 0.8 --
    ## while the proportional gradient does. The ratio cancels that level, which
    ## is what a ratio estimator is for, and the jackknife is the standard
    ## variance method for one.
    ##
    ## The absolute estimate is carried on every row as a SENSITIVITY, with its
    ## own SE and p, together with a Fieller interval for the ratio built from
    ## the joint covariance of (beta, mean). Fieller's bounded-interval condition
    ## depends on the denominator's own variance alone and ignores the strong
    ## correlation between the two halves (-0.76 to -0.85 here), so it is read as
    ## an INTERPRETABILITY flag on the percentage and never as the primary
    ## interval; on these data it is negative in every family.
    ##
    ## A family whose mean omega is at or below zero has no detectable exposure
    ## contribution and the ratio's sign is meaningless, so it reports the
    ## absolute scale. There is no threshold above zero: the old
    ## relative_scale_min_mean_z used a z that treated correlated VMRs as
    ## independent, it admitted the family it was written to exclude, and
    ## instability now shows up where it belongs -- in the ratio's own interval.
    ax <- prepare_axis(d, predictor, axis_covs)
    jt_primary <- axis_estimate_joint(ax, d$omega)
    est_absolute <- jt_primary[["beta"]]
    mean_om <- jt_primary[["mean"]]
    scale <- if (mean_om > 0) "relative_to_mean" else "raw"
    scale_reason <- if (mean_om > 0) {
        "family mean omega positive: primary is the proportional gradient"
    } else {
        "family mean omega at or below zero: no detectable exposure contribution, so the ratio has no meaning and the absolute gradient is primary"
    }
    est_primary <- if (identical(scale, "relative_to_mean")) {
        est_absolute / mean_om
    } else est_absolute
    est_descr <- axis_estimate(ax, d$y, "raw")

    ## The arm reuses `ax$ok`, so it is fitted on exactly the primary's VMRs
    ## rather than on the larger set a smaller adjustment set would admit. That
    ## is deliberate: the question is what the covariate does, not what a
    ## different denominator does.
    ax_arm <- if (arm_declared) prepare_axis(d, predictor, arm_covs,
                                            rows = ax$ok) else NULL
    jt_arm <- if (arm_declared) axis_estimate_joint(ax_arm, d$omega) else NULL
    est_arm <- if (!arm_declared) NA_real_ else if (identical(scale, "relative_to_mean")) {
        jt_arm[["beta"]] / jt_arm[["mean"]]
    } else jt_arm[["beta"]]

    ## ------------------------------------------------------------- bootstrap
    ## Donors are resampled WITHIN diagnosis so each draw keeps the observed
    ## case/control counts; inside a single-diagnosis stratum that is a plain
    ## resample. Only the outcome changes, so the axis QR is reused.
    strata_don <- as.character(d_don$diagnosis)
    ## Two columns per draw, the absolute beta and the family mean. The primary
    ## ratio series is their quotient -- identical to dividing the outcome by the
    ## draw's own mean before fitting, which is what the -b runs did -- and the
    ## two columns also give the covariance the sensitivity and Fieller need. One
    ## fit, three quantities.
    jb_p <- matrix(NA_real_, B, 2L, dimnames = list(NULL, c("beta", "mean")))
    jb_a <- matrix(NA_real_, B, 2L, dimnames = list(NULL, c("beta", "mean")))
    boot_p <- rep(NA_real_, B); boot_a <- rep(NA_real_, B)
    boot_d <- rep(NA_real_, B)
    set.seed(seed_for(opts$run_id, region, paste0("axis-boot-", ex, "-", st)))
    for (b in seq_len(B)) {
        idx <- resample_rows(strata_don)
        Xb <- X[idx, , drop = FALSE]
        if (qr(Xb)$rank < ncol(Xb)) next        # counted below, never imputed
        fmb <- fit_term_matrix(Xb, Yf[idx, , drop = FALSE], des$cols)
        om_b <- debiased_partial_ss(fmb$ss_term, fmb$df_term, fmb$sigma2, n_don)
        ## The descriptive outcome from the same draw: F = (SS/df)/sigma2.
        f_b <- (fmb$ss_term / fmb$df_term) / fmb$sigma2
        p_b <- stats::pf(f_b, fmb$df_term, fmb$df_resid, lower.tail = FALSE)
        boot_d[b] <- axis_estimate(ax, -log10(pmax(p_b, eps)), "raw")
        jb_p[b, ] <- axis_estimate_joint(ax, om_b)
        if (arm_declared) jb_a[b, ] <- axis_estimate_joint(ax_arm, om_b)
        ## The ratio series needs a positive denominator in the draw; such a
        ## draw is counted in n_bootstrap_failed, never imputed. The absolute
        ## series above keeps every draw, since it has no denominator.
        if (identical(scale, "relative_to_mean")) {
            if (jb_p[b, "mean"] > 0) boot_p[b] <- jb_p[b, "beta"] / jb_p[b, "mean"]
            if (arm_declared && jb_a[b, "mean"] > 0) {
                boot_a[b] <- jb_a[b, "beta"] / jb_a[b, "mean"]
            }
        } else {
            boot_p[b] <- jb_p[b, "beta"]
            if (arm_declared) boot_a[b] <- jb_a[b, "beta"]
        }
    }
    ## Primary: the jackknife of the SAME functional the point estimate is, so a
    ## deleted block recomputes the ratio on the rows it keeps.
    se_jk_p <- block_jackknife_se(ax, d$omega, d$chrom, scale)
    se_jk_d <- block_jackknife_se(ax, d$y, d$chrom, "raw")
    inf_p <- combined_inference(est_primary, boot_p, se_jk_p, alpha_ci)
    inf_d <- combined_inference(est_descr, boot_d, se_jk_d, alpha_ci)
    ## Sensitivity: the absolute gradient and the Fieller flag, from the joint
    ## covariance of the two halves.
    cov_jk_p <- block_jackknife_cov(ax, d$omega, d$chrom)
    V_p <- combined_cov(jb_p, cov_jk_p)
    inf_abs <- combined_inference(est_absolute, jb_p[, "beta"],
                                  sqrt(cov_jk_p["beta", "beta"]), alpha_ci)
    rel_p <- fieller_ratio_ci(est_absolute, mean_om, V_p, alpha_ci)
    ## How far the donor bootstrap inflates the denominator. Resampling with
    ## replacement duplicates donors, which inflates an exposure-explained sum of
    ## squares; measured at 2.1x, 5.0x and 38x on the -b runs. It is on every row
    ## because it is the reason the ratio's bootstrap variance may be optimistic,
    ## and that caveat should be auditable from the table rather than only in
    ## prose. AGENTS.md 7.9 already bars a percentile interval from this
    ## bootstrap; only its variance is used.
    ## Only meaningful where the ratio is the primary: with a mean at or below
    ## zero the quotient of two means has no interpretation as an inflation.
    boot_mean_inflation <- if (identical(scale, "relative_to_mean")) {
        mean(jb_p[, "mean"], na.rm = TRUE) / mean_om
    } else NA_real_
    cov_jk_a <- if (arm_declared) block_jackknife_cov(ax_arm, d$omega, d$chrom) else NULL
    inf_a <- if (arm_declared) {
        combined_inference(est_arm, boot_a,
                           block_jackknife_se(ax_arm, d$omega, d$chrom, scale),
                           alpha_ci)
    } else NULL
    rel_a <- if (arm_declared) {
        fieller_ratio_ci(jt_arm[["beta"]], jt_arm[["mean"]],
                         combined_cov(jb_a, cov_jk_a), alpha_ci)
    } else NULL

    ## The OLS p-value this stage used to report, kept for the comparison only.
    ols <- summary(stats::lm(stats::as.formula(paste("omega ~", predictor, "+",
                                                    covar_rhs)), data = d))
    ols_d <- summary(stats::lm(stats::as.formula(paste("y ~", predictor, "+",
                                                      covar_rhs)), data = d))

    ## Threshold-based secondaries, unchanged.
    w <- if (n_sig >= 10L && n_sig <= n_axis - 10L) {
        stats::wilcox.test(d[[predictor]][d$significant],
                           d[[predictor]][!d$significant], exact = FALSE)
    } else NULL
    g <- if (n_sig >= 10L && n_sig <= n_axis - 10L) {
        f3 <- stats::as.formula(paste("significant ~", predictor, "+", covar_rhs))
        fit <- try(stats::glm(f3, data = d, family = stats::binomial()), silent = TRUE)
        if (inherits(fit, "try-error")) NULL else fit
    } else NULL
    g_co <- if (!is.null(g)) summary(g)$coefficients[predictor, ] else NULL

    cbind(base, data.table(
        status = "ok",
        n_vmrs_in_axis = n_axis,
        n_donors_refit = n_don,
        primary_outcome = "debiased_partial_ss",
        primary_outcome_meaning = paste(
            "proportional change in the debiased exposure-explained sum of",
            "squares per donor, per SD of score, as a fraction of this family's",
            "mean; multiply by 100 for a percentage. Families whose mean is at",
            "or below zero report the absolute gradient instead, and say so in",
            "primary_scale_reason"),
        primary_scale = scale,
        primary_scale_reason = scale_reason,
        mean_omega = mean_om,
        mean_omega_se = sqrt(V_p["mean", "mean"]),
        mean_omega_z = rel_p$den_z,
        frac_omega_positive = mean(d$omega > 0),
        primary_beta = est_primary,
        primary_se = inf_p$se,
        primary_se_bootstrap = inf_p$se_bootstrap,
        primary_se_jackknife = inf_p$se_jackknife,
        primary_ci_lower = inf_p$ci_lower,
        primary_ci_upper = inf_p$ci_upper,
        primary_p = inf_p$p,
        primary_p_ols_vmrs_independent = unname(ols$coefficients[predictor, 4]),
        n_bootstrap = B,
        n_bootstrap_used = inf_p$n_bootstrap_used,
        n_bootstrap_failed = B - inf_p$n_bootstrap_used,
        primary_model = paste0("omega ~ ", predictor, " + ", covar_rhs),
        inference = "donor_bootstrap_var_plus_chromosome_jackknife_var",
        ## ------------------------------------------------- absolute sensitivity
        ## The same gradient in the outcome's own units, with the level of
        ## exposure-explained variance left in rather than divided out. It is a
        ## SENSITIVITY: it answers a different question -- does the absolute
        ## gradient replicate across chromosomes -- and on these data the answer
        ## is no, because chromosomes differ about five-fold in that level. Its
        ## block-jackknife SE runs 5-10x its bootstrap SE for that reason. It
        ## carries no FDR family and can neither promote nor demote the primary.
        absolute_beta = est_absolute,
        absolute_se = inf_abs$se,
        absolute_se_bootstrap = inf_abs$se_bootstrap,
        absolute_se_jackknife = inf_abs$se_jackknife,
        absolute_ci_lower = inf_abs$ci_lower,
        absolute_ci_upper = inf_abs$ci_upper,
        absolute_p = inf_abs$p,
        absolute_role = "sensitivity_not_gating",
        ## Fieller's interval for the ratio, built from the joint covariance of
        ## (absolute beta, mean omega). Informational only. Its bounded-interval
        ## condition depends on the denominator's own variance and ignores the
        ## strong negative correlation between the two halves, so it is
        ## conservative here to the point of being uninformative: it returns
        ## not-estimable in every family, including the ones the primary
        ## resolves. Recorded so that conservatism is on the record rather than
        ## asserted, and never used as a gate.
        fieller_ci_lower = rel_p$ci_lower,
        fieller_ci_upper = rel_p$ci_upper,
        fieller_bounded = rel_p$estimable,
        fieller_unbounded_reason = rel_p$reason,
        fieller_role = "informational_not_gating",
        ## Donor-bootstrap inflation of the denominator: resampling with
        ## replacement duplicates donors, which inflates an exposure-explained
        ## sum of squares. Above 1 the ratio's bootstrap variance is likely
        ## optimistic, so this is the audit trail for that caveat.
        bootstrap_mean_omega_inflation = boot_mean_inflation,
        ## Descriptive, and never a decision. Its expectation depends on the SE,
        ## so it is biased toward the hypothesis by construction.
        ## Non-gating arm: the primary model minus methylation_variance.
        ## `arm_attenuation` is the fraction of the primary coefficient that the
        ## covariate carries; it is NOT a test, and neither sign is a failure.
        arm_no_methylation_variance = arm_declared,
        arm_covariates = if (arm_declared) paste(arm_covs, collapse = ",") else NA_character_,
        arm_beta = if (arm_declared) est_arm else NA_real_,
        arm_se = if (arm_declared) inf_a$se else NA_real_,
        arm_p = if (arm_declared) inf_a$p else NA_real_,
        arm_ci_lower = if (arm_declared) inf_a$ci_lower else NA_real_,
        arm_ci_upper = if (arm_declared) inf_a$ci_upper else NA_real_,
        arm_fieller_bounded = if (arm_declared) rel_a$estimable else NA,
        arm_attenuation = if (arm_declared && is.finite(est_primary) &&
                              abs(est_primary) > 0) 1 - est_arm / est_primary else NA_real_,
        arm_role = "non_gating_sensitivity",
        neglog10p_beta = est_descr,
        neglog10p_se = inf_d$se,
        neglog10p_p = inf_d$p,
        neglog10p_p_ols_vmrs_independent = unname(ols_d$coefficients[predictor, 4]),
        mechanically_biased_toward_hypothesis = TRUE,
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

## Correction is WITHIN STRATUM here too, matching Stage A. The pooled and
## within-case analyses run on different donor sets and answer different
## questions, so pooling them into one family would mix them exactly the way
## AGENTS.md 10.3 forbids -- and the family would then depend on how many
## exposures happened to clear the gate in the *other* stratum.
##
## The three model types are also separate families: the primary is
## threshold-free and the two grouped tests are secondary, so a secondary test
## must not borrow significance from the primary or vice versa.
out[status == "ok",
    primary_fdr := stats::p.adjust(primary_p, method = "BH"), by = stratum]
out[status == "ok",
    neglog10p_fdr := stats::p.adjust(neglog10p_p, method = "BH"), by = stratum]
## Its own family: a non-gating sensitivity must not borrow significance from the
## primary, nor lend it.
out[status == "ok" & !is.na(arm_p),
    arm_fdr := stats::p.adjust(arm_p, method = "BH"), by = stratum]
out[status == "ok" & !is.na(wilcoxon_p),
    wilcoxon_fdr := stats::p.adjust(wilcoxon_p, method = "BH"), by = stratum]
out[status == "ok" & !is.na(logistic_p),
    logistic_fdr := stats::p.adjust(logistic_p, method = "BH"), by = stratum]

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
    axis_covariates = paste(axis_covs, collapse = ","),
    axis_primary_outcome = "debiased_partial_ss",
    axis_primary_scale = "proportional_ratio_of_family_mean",
    axis_absolute_sensitivity = "absolute_beta_with_fieller_flag",
    n_absolute_p_below_alpha = as.character(
        nrow(out[status == "ok" & absolute_p < alpha_ci])),
    n_fieller_bounded = as.character(nrow(out[status == "ok" & fieller_bounded == TRUE])),
    max_bootstrap_mean_omega_inflation = {
        v <- out$bootstrap_mean_omega_inflation
        v <- v[is.finite(v)]
        if (length(v)) as.character(signif(max(v), 4)) else NA_character_
    },
    axis_non_gating_arms = if (arm_declared) "no_methylation_variance" else "none",
    axis_inference = "donor_bootstrap_var_plus_chromosome_jackknife_var",
    n_bootstrap = as.character(B),
    n_axis_donors_refit = paste(sort(unique(out$n_donors_refit)), collapse = ",")
))
print(out[, .(exposure, stratum, primary_scale, primary_beta, primary_p,
              primary_fdr, mean_omega, absolute_beta, absolute_p,
              bootstrap_mean_omega_inflation, arm_beta, arm_p, neglog10p_p)])
