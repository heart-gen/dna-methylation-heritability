## Fold-internal Haseman-Elston screen, calibrated by permutation.
##
## AGENTS.md 7.3: "Do not use the current ordinary-OLS Haseman-Elston p-value as
## the production screen without donor-robust or permutation calibration." The
## p-value that `haseman_elston()` returns is exactly that forbidden quantity:
## it treats the n*(n-1)/2 pairwise products as independent observations, which
## they are not, so it is anticonservative by a large and locus-dependent
## factor. This file supplies the permutation calibration instead.
##
## Everything here is computed on OUTER-TRAINING donors only.
##
## Why this is not simply a loop over `haseman_elston()`:
##
##   The HE slope regresses the pairwise product v_ij = r_i r_j on the pairwise
##   relatedness u_ij = GRM_ij, over the upper triangle, with an intercept:
##
##       slope = sum_upper (u - ubar) v / sum_upper (u - ubar)^2
##
##   Under permutation only r changes; the GRM, and therefore the entire
##   denominator, is fixed. Writing C for the symmetric matrix with
##   C_ij = u_ij - ubar off the diagonal and 0 on it,
##
##       sum_upper (u - ubar) v = 0.5 * r' C r
##
##   so every permutation replicate is a quadratic form in r. Stacking B
##   permuted residual vectors as the columns of R turns the whole null
##   distribution into ONE matrix product, colSums(R * (C %*% R)). At n = 153
##   and B = 1000 that is roughly 20 MFLOP through BLAS -- milliseconds --
##   whereas B calls to `haseman_elston()` would each rebuild an n x n GRM and
##   an n(n-1)/2-row design matrix.
##
##   `he_slope_equals_haseman_elston` in tests/ asserts the identity on random
##   data, so this stays a fast path to the SAME statistic rather than a second,
##   silently diverging estimator.

#' Precompute the genotype-only half of the screen for one outer-training fold.
#'
#' Returns NULL when the locus has too little genotype variation to screen,
#' which the caller must treat as a screen FAILURE (null prediction), not an
#' error.
he_screen_prepare <- function(g_train) {
    g <- as.matrix(g_train)
    storage.mode(g) <- "double"
    variances <- apply(g, 2L, stats::var)
    keep <- is.finite(variances) & variances > 1e-8
    g <- g[, keep, drop = FALSE]
    if (ncol(g) < 2L) return(NULL)

    standardized <- scale(g)
    grm <- tcrossprod(standardized) / ncol(standardized)

    upper <- upper.tri(grm)
    u <- grm[upper]
    ubar <- mean(u)
    denom <- sum((u - ubar)^2)
    if (!is.finite(denom) || denom <= 0) return(NULL)

    ## C: relatedness centered off-diagonal, zeroed on the diagonal, so that
    ## 0.5 * r' C r reproduces the upper-triangle sum exactly.
    cmat <- grm - ubar
    diag(cmat) <- 0
    list(cmat = cmat, denom = denom, n_snps = ncol(g), n = nrow(g))
}

#' HE slope for one residual vector, given a prepared fold.
he_slope <- function(prep, r) {
    r <- as.numeric(scale(r))
    0.5 * drop(crossprod(r, prep$cmat %*% r)) / prep$denom
}

#' Permutation-calibrated fold-internal screen.
#'
#' @param prep from he_screen_prepare()
#' @param y_train training phenotype, ALREADY residualized on covariates in
#'   this fold. Permuting a covariate-residualized phenotype is what makes the
#'   exchangeability assumption defensible here.
#' @param n_perm permutation replicates (config screen.n_permutations)
#' @param alpha threshold (config screen.alpha)
#' @param seed deterministic, from seed_for()
#'
#' One-sided on purpose: the alternative is positive local genetic control, and
#' a strongly negative HE slope is noise, not evidence for prediction.
#'
#' p = (1 + #{perm >= observed}) / (1 + n_perm) -- the add-one form, so p is
#' never 0 and the screen cannot claim more resolution than n_perm supports.
he_permutation_screen <- function(prep, y_train, n_perm, alpha, seed) {
    if (is.null(prep)) {
        return(list(pass = FALSE, p = NA_real_, statistic = NA_real_,
                    n_perm = 0L, reason = "insufficient_genotype_variation"))
    }
    observed <- he_slope(prep, y_train)
    if (!is.finite(observed)) {
        return(list(pass = FALSE, p = NA_real_, statistic = NA_real_,
                    n_perm = 0L, reason = "nonfinite_he_statistic"))
    }

    n <- length(y_train)
    r0 <- as.numeric(scale(y_train))
    set.seed(seed)
    perm <- matrix(r0[replicate(n_perm, sample.int(n))], nrow = n)
    ## Columns are already standardized: permutation only reorders r0.
    stats_null <- 0.5 * colSums(perm * (prep$cmat %*% perm)) / prep$denom

    finite_null <- stats_null[is.finite(stats_null)]
    if (length(finite_null) == 0L) {
        return(list(pass = FALSE, p = NA_real_, statistic = observed,
                    n_perm = 0L, reason = "no_finite_permutation_statistic"))
    }
    p <- (1 + sum(finite_null >= observed)) / (1 + length(finite_null))
    list(pass = p <= alpha, p = p, statistic = observed,
         n_perm = length(finite_null), reason = NA_character_)
}
