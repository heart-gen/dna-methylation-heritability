#### Shared axis-model estimation and inference ####
##
## An "axis model" regresses a per-VMR outcome on Module 02's
## local_snp_contribution_score_z plus technical covariates, over the VMRs of one
## cohort x region cell. Modules 09b and 10 both do this, and both hit the same
## two problems, so the machinery lives here rather than in either module
## (AGENTS.md 5.3).
##
## PROBLEM 1 -- the outcome scale. The per-VMR model that produces the outcome
## (an age slope in 09b, an exposure effect in 10) contains no SNP term, so a VMR
## with strong local SNP control carries its local genetic variance in that
## model's residual. Its SE is therefore larger, and ANY outcome whose
## expectation depends on the SE -- |beta_hat|, its rank, |t|, -log10 p --
## couples to the score mechanically. The debiased quantities below have
## expectations that do not involve the SE:
##
##   debiased_sq_effect(beta, se)               = beta_hat^2 - SE^2
##   debiased_partial_ss(ss, df, sigma2, n)     = (SS_exposure - df*sigma2)/n
##
## The second is the first generalized to a multi-df term. SS_exposure is the
## extra sum of squares of the term against the covariate-only model, which
## equals beta_E' M^-1 beta_E with M the term's block of (X'X)^-1, so
## E[SS] = beta_E' M^-1 beta_E + df*sigma^2 and subtracting df*sigma2_hat leaves
## an unbiased estimate of the term's contribution. For a 1-df term the two
## agree up to a factor that is constant within a family, because the design is
## shared by every VMR in it. Both are zero in expectation under the null and
## both take negative values, which are retained, never truncated.
##
## PROBLEM 2 -- the variance. An axis model has 9,000-11,000 rows, but the VMRs
## are neither independent of each other (nearby loci covary) nor independent
## across the donor sample (every VMR is estimated in the same ~120-150 donors).
## An OLS SE over the VMR rows treats both as independent and understates the
## uncertainty; in Module 10 it ran a median 1.28x small against the block
## jackknife alone. Inference is therefore
##
##   SE^2 = donor-bootstrap variance + delete-one-chromosome jackknife variance
##
## and either half alone under-covers. The bootstrap supplies a VARIANCE only:
## never form a percentile interval from the bootstrap of a debiased outcome,
## whose bootstrap distribution sits about SE^2 away from the estimate by
## construction.
##
## Provenance: extracted from 09b_aging_application/_h/age_functions.R on
## 2026-09-19, arithmetic unchanged, when Module 10 needed the same three
## functions. 09b sources this file and keeps only its age-specific parts.

## ---------------------------------------------------------------- debiasing

debiased_sq_effect <- function(beta, se) beta^2 - se^2

debiased_partial_ss <- function(ss_term, df_term, sigma2, n) {
    (ss_term - df_term * sigma2) / n
}

## ------------------------------------------------- per-VMR matrix estimation

## One QR per design applied to every VMR at once. `col` is the index of the
## predictor of interest in X, and the pivot is honoured because qr() may
## reorder columns. Returns beta, SE, t and p for that one column.
fit_scalar_matrix <- function(X, Y, col) {
    q <- qr(X)
    beta_all <- qr.coef(q, Y)
    resid <- qr.resid(q, Y)
    df <- nrow(X) - q$rank
    sigma2 <- colSums(resid^2) / df
    xtx_inv <- chol2inv(qr.R(q))
    a <- which(q$pivot == col)
    se <- sqrt(sigma2 * xtx_inv[a, a])
    beta <- if (is.matrix(beta_all)) beta_all[col, ] else beta_all[col]
    t <- beta / se
    list(beta = as.numeric(beta), se = as.numeric(se), t = as.numeric(t),
         p = as.numeric(2 * stats::pt(-abs(t), df)), sigma2 = as.numeric(sigma2),
         df = df)
}

## The multi-df analogue: the extra sum of squares of the columns `cols` of X,
## for every VMR at once, plus the full model's sigma2. Computed as
## RSS(covariates only) - RSS(full), which is what anova() reports and needs no
## explicit (X'X)^-1 block.
fit_term_matrix <- function(X, Y, cols) {
    q_full <- qr(X)
    rss_full <- colSums(qr.resid(q_full, Y)^2)
    df_full <- nrow(X) - q_full$rank
    X0 <- X[, setdiff(seq_len(ncol(X)), cols), drop = FALSE]
    q0 <- qr(X0)
    rss_null <- colSums(qr.resid(q0, Y)^2)
    list(ss_term = as.numeric(rss_null - rss_full),
         sigma2 = as.numeric(rss_full / df_full),
         df_term = q_full$rank - q0$rank,
         df_resid = df_full)
}

## ------------------------------------------------------------- donor sampling

## Resample donors with replacement WITHIN strata, so every draw keeps the
## observed group counts and the design stays estimable.
resample_rows <- function(strata) {
    idx <- integer(0)
    for (s in unique(strata)) {
        pool <- which(strata == s)
        idx <- c(idx, pool[sample.int(length(pool), length(pool), replace = TRUE)])
    }
    idx
}

## ------------------------------------------------------------- the axis model

## Axis design for a fixed set of VMR rows. Only the outcome changes between
## bootstrap draws, so the QR is taken once and every refit is a single qr.coef.
prepare_axis <- function(dt, predictor, covariates, rows = NULL) {
    cols <- c(predictor, covariates)
    ok <- stats::complete.cases(dt[, cols, with = FALSE])
    if (!is.null(rows)) ok <- ok & rows
    Z <- cbind(`(Intercept)` = 1, as.matrix(dt[ok, cols, with = FALSE]))
    storage.mode(Z) <- "double"
    q <- qr(Z)
    if (q$rank < ncol(Z)) {
        stop("Axis design is rank-deficient for covariates: ",
             paste(covariates, collapse = ","))
    }
    list(ok = ok, Z = Z, qr = q, j = which(colnames(Z) == predictor),
         predictor = predictor, covariates = covariates)
}

## The axis estimate for one outcome on a subset of the axis rows. A
## "relative_to_mean" outcome is divided by its mean over the SAME rows, so the
## coefficient reads as the proportional change in the mean outcome per SD of
## score -- unit-free, and so comparable between cells whose absolute scales
## differ. `drop` removes rows (the jackknife).
axis_estimate <- function(ax, y, scale = c("raw", "relative_to_mean"),
                          drop = NULL) {
    scale <- match.arg(scale)
    yy <- y[ax$ok]
    Z <- ax$Z
    if (!is.null(drop)) {
        keep <- !drop[ax$ok]
        yy <- yy[keep]; Z <- Z[keep, , drop = FALSE]
        q <- qr(Z)
    } else {
        q <- ax$qr
    }
    if (scale == "relative_to_mean") yy <- yy / mean(yy)
    unname(qr.coef(q, yy)[ax$j])
}

## Delete-one-chromosome weighted block jackknife (Busing, Meijer & van der
## Leeden 1999) -- the construction Modules 06 and 08 use. Weighted because
## chromosomes differ several-fold in VMR count. This is the VMR-level half of
## the variance: the VMRs are a sample of loci, and nearby loci are not
## independent, so the block is the chromosome.
block_jackknife_se <- function(ax, y, chrom, scale) {
    blocks <- sort(unique(chrom[ax$ok]))
    if (length(blocks) < 2L) return(NA_real_)
    full <- axis_estimate(ax, y, scale)
    n_tot <- sum(ax$ok)
    hj <- n_tot / vapply(blocks, function(b) sum(chrom[ax$ok] == b), numeric(1))
    theta <- vapply(blocks, function(b) axis_estimate(ax, y, scale,
                                                      drop = chrom == b), numeric(1))
    pseudo <- hj * full - (hj - 1) * theta
    sqrt(sum((pseudo - mean(pseudo))^2 / (hj - 1)) / length(blocks))
}

## Combined inference. Donor-bootstrap variance carries the donor-level half
## (shared donors correlate every VMR's estimate); the block jackknife carries
## the VMR-level half. Adding them double-counts the part of the per-VMR noise
## both see, which errs conservative; in 09b's simulations the donor bootstrap
## alone rejected 57% of true nulls, the sum 0-3%.
combined_inference <- function(estimate, boot, se_jk, alpha = 0.05) {
    boot <- boot[is.finite(boot)]
    se <- sqrt(stats::var(boot) + se_jk^2)
    z <- stats::qnorm(1 - alpha / 2)
    list(se = se, se_bootstrap = stats::sd(boot), se_jackknife = se_jk,
         z = estimate / se, p = 2 * stats::pnorm(-abs(estimate / se)),
         ci_lower = estimate - z * se, ci_upper = estimate + z * se,
         n_bootstrap_used = length(boot))
}

## ---------------------------------------------------------------------------
## RATIO-AWARE RELATIVE EFFECT (PI, 2026-09-20)
##
## The absolute coefficient is in the outcome's own units (for Module 10, a
## debiased sum of squares per donor, order 1e-5), which is not a quantity
## anyone can read. The interpretable form is the fraction of the family mean,
## R = beta / mean(y). Dividing the OUTCOME by its mean before fitting -- the
## old `scale = "relative_to_mean"` path -- makes that ratio the estimand, and
## then every bootstrap draw and every deleted chromosome carries a different
## denominator. Where the mean sits near zero the coefficient explodes:
## env-AA-dlpfc-20260919-{a,b} reported -2.41 (a -241% gradient) on a mean of
## 4.9e-06, with a jackknife SE of 1.49.
##
## So estimate on the absolute scale and treat the ratio as a derived quantity
## with its own uncertainty. Fieller's theorem inverts the test of
## beta - R * mean = 0 over R, using the joint covariance of the two estimates,
## and it returns no bounded interval exactly when the denominator is not itself
## separated from zero at the same level. That replaces a hand-set stability
## threshold with the confidence level already in use: the data decide whether a
## percentage is estimable.

## Both halves of the ratio, from one fit on one row set.
axis_estimate_joint <- function(ax, y, drop = NULL) {
    yy <- y[ax$ok]
    Z <- ax$Z
    if (!is.null(drop)) {
        keep <- !drop[ax$ok]
        yy <- yy[keep]; Z <- Z[keep, , drop = FALSE]
        q <- qr(Z)
    } else {
        q <- ax$qr
    }
    c(beta = unname(qr.coef(q, yy)[ax$j]), mean = mean(yy))
}

## The jackknife of both halves together, so the COVARIANCE is carried and not
## just the two variances. Same weighted construction as block_jackknife_se();
## the scalar function is left alone because 09b's accepted runs call it.
block_jackknife_cov <- function(ax, y, chrom) {
    blocks <- sort(unique(chrom[ax$ok]))
    if (length(blocks) < 2L) return(matrix(NA_real_, 2L, 2L))
    full <- axis_estimate_joint(ax, y)
    n_tot <- sum(ax$ok)
    hj <- n_tot / vapply(blocks, function(b) sum(chrom[ax$ok] == b), numeric(1))
    theta <- vapply(blocks, function(b) axis_estimate_joint(ax, y, drop = chrom == b),
                    numeric(2))                      # 2 x n_blocks
    pseudo <- outer(full, hj) - sweep(theta, 2, hj - 1, `*`)
    centred <- pseudo - rowMeans(pseudo)
    v <- tcrossprod(sweep(centred, 2, sqrt(hj - 1), `/`)) / length(blocks)
    dimnames(v) <- list(c("beta", "mean"), c("beta", "mean"))
    v
}

## Donor-bootstrap covariance plus block-jackknife covariance, the matrix form
## of combined_inference()'s variance sum and justified by the same argument.
combined_cov <- function(boot_mat, cov_jk) {
    ok <- stats::complete.cases(boot_mat)
    v <- stats::cov(boot_mat[ok, , drop = FALSE]) + cov_jk
    dimnames(v) <- list(c("beta", "mean"), c("beta", "mean"))
    v
}

## Fieller (1954) interval for num/den. Solving
##   (num - R*den)^2 = z^2 * (v11 - 2R*v12 + R^2*v22)
## gives A*R^2 + B*R + C = 0 with A = den^2 - z^2*v22. A <= 0 means the
## denominator is not distinguishable from zero at this level, so the solution
## set is unbounded (the complement of an interval) and no percentage is
## reportable. `estimable` is that verdict; it is not a significance test of the
## ratio, whose null R = 0 is identical to num = 0 and is already tested there.
fieller_ratio_ci <- function(num, den, V, alpha = 0.05) {
    z <- stats::qnorm(1 - alpha / 2)
    ratio <- num / den
    v11 <- V[1, 1]; v12 <- V[1, 2]; v22 <- V[2, 2]
    A <- den^2 - z^2 * v22
    Bq <- -2 * (num * den - z^2 * v12)
    C <- num^2 - z^2 * v11
    disc <- Bq^2 - 4 * A * C
    bad <- function(why) list(ratio = ratio, ci_lower = NA_real_,
                              ci_upper = NA_real_, estimable = FALSE,
                              reason = why,
                              den_z = den / sqrt(v22))
    if (!all(is.finite(c(num, den, v11, v12, v22)))) return(bad("nonfinite_inputs"))
    if (den <= 0) return(bad("family_mean_at_or_below_zero"))
    if (A <= 0) return(bad("family_mean_not_separated_from_zero_at_this_level"))
    if (disc < 0) return(bad("fieller_solution_set_empty"))
    lo <- (-Bq - sqrt(disc)) / (2 * A)
    hi <- (-Bq + sqrt(disc)) / (2 * A)
    list(ratio = ratio, ci_lower = min(lo, hi), ci_upper = max(lo, hi),
         estimable = TRUE, reason = NA_character_, den_z = den / sqrt(v22))
}
