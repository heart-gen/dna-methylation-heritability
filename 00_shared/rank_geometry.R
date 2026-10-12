#### Locus-geometry decomposition of the frozen joint estimate ####
##
## Not sourced by load.R. Scripts that need it source it explicitly after
## load.R:
##
##     source(file.path(repo_root(), "00_shared", "rank_geometry.R"))
##
## WHY THIS EXISTS
## The frozen joint model of 02_local_genetic_variance turns four estimator
## outputs into the unbounded estimate that Stage 04 ranks. It also takes four
## locus-geometry terms -- log n, log num_snps, log p_eff and ld_metric -- and
## enters them LINEARLY (joint_pve_functions.R::make_joint_matrix). Within one
## cell n is constant, so three geometry terms move a VMR's estimate up or down
## regardless of its methylation data. If SNP count, effective rank or LD differ
## between chromatin states, those terms alone could produce a chromatin
## association with the rank. This file isolates them so that question can be
## answered on the sealed estimates.
##
## HOW
## The model is linear in (a) monotone hinge expansions of bslmm_pve, he_h2,
## rho2_oof and r2_oof at the fixed knots and clipping bounds of
## config/joint-pve-20260820.tsv, and (b) the geometry terms. Regressing the
## sealed pve_cis_joint_unbounded on that same basis by OLS therefore recovers
## the model exactly, without the fitted object. Two checks make that claim
## testable rather than assumed, and the run stops if either fails:
##
##   * the recovery R^2 is 1 to within `tol` in every cell;
##   * the recovered geometry coefficients agree across cells, because one
##     frozen model scored them all.
##
## The geometry-free estimate is the recovered fit minus its geometry part. It
## is ranked with Stage 04's own rule (midrank percentile, then standardized),
## and the same rule applied to the sealed estimate must reproduce the sealed
## local_snp_contribution_score_z, which checks the rule as well.
##
## WHAT THIS IS NOT
## The geometry-free rank is a sensitivity predictor. It is not a replacement
## score, carries no calibration of its own, and must not be reported as a
## locus-level quantity.

JOINT_SIGNALS <- c("bslmm_pve", "he_h2", "rho2_oof", "r2_oof")
GEOMETRY_TERMS <- c(log_num_snps = "num_snps", log_p_eff = "p_eff",
                    ld_metric = "ld_metric")

#' Read the frozen joint model's settings table (setting/value TSV).
read_joint_settings_table <- function(path) {
    tab <- data.table::fread(path, colClasses = "character")
    if (!all(c("setting", "value") %in% names(tab))) {
        stop("Joint settings table lacks setting/value columns: ", path)
    }
    stats::setNames(as.list(tab$value), tab$setting)
}

.split_num <- function(x) as.numeric(trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1L]]))

#' The hinge basis of the four estimator outputs, exactly as the frozen model
#' builds it: clip each output to its bounds, then max(x - knot, 0) per knot.
joint_hinge_basis <- function(dt, settings) {
    cols <- list()
    for (f in JOINT_SIGNALS) {
        ## Settings keys use the output's prefix: bslmm_clip, he_knots, ...
        prefix <- sub("_.*$", "", f)
        clip <- .split_num(settings[[paste0(prefix, "_clip")]])
        knots <- .split_num(settings[[paste0(prefix, "_knots")]])
        if (length(clip) != 2L || !all(is.finite(clip)) || !length(knots)) {
            stop("Joint settings lack clip/knots for ", f)
        }
        x <- as.numeric(dt[[f]])
        if (any(!is.finite(x))) stop("Nonfinite estimator output in ", f)
        x <- pmin(clip[[2L]], pmax(clip[[1L]], x))
        for (k in knots) cols[[sprintf("%s__k%s", f, format(k))]] <- pmax(x - k, 0)
    }
    do.call(cbind, cols)
}

#' Geometry terms on the scale the frozen model uses (log for counts).
joint_geometry_terms <- function(dt) {
    g <- cbind(log_num_snps = log(as.numeric(dt$num_snps)),
               log_p_eff = log(as.numeric(dt$p_eff)),
               ld_metric = as.numeric(dt$ld_metric))
    if (any(!is.finite(g))) stop("Nonfinite geometry term")
    g
}

#' Stage 04's within-cell rule: midrank percentile, then standardized.
within_cell_rank_z <- function(x) {
    s <- (rank(x, ties.method = "average") - 0.5) / length(x)
    as.numeric(scale(s))
}

#' Decompose one cell's eligible VMRs.
#'
#' @param dt eligible rows of a sealed Module 02 per-VMR table
#' @return list(per_vmr, coefficients, summary)
decompose_rank_geometry <- function(dt, settings, cell, region, tol = 1e-8) {
    need <- c("vmr_id", JOINT_SIGNALS, unname(GEOMETRY_TERMS),
              "pve_cis_joint_unbounded", "local_snp_contribution_score_z")
    miss <- setdiff(need, names(dt))
    if (length(miss)) stop(cell, " ", region, ": table lacks ", paste(miss, collapse = ", "))
    y <- as.numeric(dt$pve_cis_joint_unbounded)
    H <- joint_hinge_basis(dt, settings)
    G <- joint_geometry_terms(dt)
    fit <- stats::lm.fit(cbind(`(Intercept)` = 1, H, G), y)
    b <- fit$coefficients
    b[is.na(b)] <- 0
    fitted <- y - fit$residuals
    r2 <- 1 - sum(fit$residuals^2) / sum((y - mean(y))^2)
    if (!is.finite(r2) || r2 < 1 - tol) {
        stop(cell, " ", region, ": the hinge + geometry basis recovers the sealed ",
             "estimate at R^2 = ", format(r2, digits = 12), ", not 1. The basis no ",
             "longer matches the frozen model; nothing below would be about it.")
    }
    geom <- drop(G %*% b[colnames(G)])
    geom <- geom - mean(geom)
    signal <- fitted - geom

    z_check <- within_cell_rank_z(y)
    z_err <- max(abs(z_check - as.numeric(dt$local_snp_contribution_score_z)))
    if (!is.finite(z_err) || z_err > 1e-8) {
        stop(cell, " ", region, ": Stage 04's rank rule does not reproduce the sealed ",
             "local_snp_contribution_score_z (max abs diff ", format(z_err, digits = 3), ")")
    }
    z_free <- within_cell_rank_z(signal)
    sp <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman"))

    per_vmr <- data.table::data.table(
        cell = cell, region = region, vmr_id = dt$vmr_id,
        local_snp_contribution_score_z = as.numeric(dt$local_snp_contribution_score_z),
        geometry_component = geom,
        geometry_free_estimate = signal,
        geometry_free_rank_z = z_free)
    coefficients <- data.table::data.table(
        cell = cell, region = region, term = colnames(G),
        coefficient = unname(b[colnames(G)]))
    summary <- data.table::data.table(
        cell = cell, region = region, n_eligible = nrow(dt),
        recovery_r2 = r2,
        rank_rule_max_abs_diff = z_err,
        var_share_geometry = stats::var(geom) / stats::var(y),
        spearman_rank_vs_geometry_free_rank = sp(y, signal),
        spearman_rank_vs_geometry_component = sp(y, geom),
        spearman_rank_vs_num_snps = sp(y, dt$num_snps),
        spearman_rank_vs_p_eff = sp(y, dt$p_eff),
        spearman_rank_vs_ld_metric = sp(y, dt$ld_metric),
        spearman_geometry_free_rank_vs_num_snps = sp(signal, dt$num_snps),
        spearman_geometry_free_rank_vs_p_eff = sp(signal, dt$p_eff),
        spearman_geometry_free_rank_vs_ld_metric = sp(signal, dt$ld_metric))
    list(per_vmr = per_vmr, coefficients = coefficients, summary = summary)
}

#' Assert one frozen model: the recovered geometry coefficients agree across
#' cells to `tol` (relative to the largest coefficient).
assert_frozen_geometry <- function(coefficients, tol = 1e-6) {
    spread <- coefficients[, .(spread = max(coefficient) - min(coefficient),
                               scale = max(abs(coefficient))), by = term]
    bad <- spread[spread > tol * pmax(scale, 1)]
    if (nrow(bad)) {
        print(bad)
        stop("Recovered geometry coefficients differ between cells; the cells were ",
             "not scored by one frozen model, or the basis is wrong.")
    }
    invisible(spread)
}
