#!/usr/bin/env Rscript

## Stage 16: the reliability ceiling for the donor-group concordance result.
##
## WHY THIS EXISTS
## Stage 14 asks whether the ORDERING of loci agrees between donor groups. A
## rank correlation from that table is uninterpretable on its own, because two
## different things depress it and they have opposite meanings:
##
##   1. the donor groups genuinely differ in local genetic control, or
##   2. neither cell can reproduce its OWN ordering, so there is nothing for
##      the two to agree about.
##
## Spearman 0.55 is strong evidence of shared architecture if each cell can only
## reproduce itself at 0.6, and evidence of nothing if each reproduces itself at
## 0.95. The ceiling is what separates those readings, and without it the
## concordance number should not be reported at all.
##
## WHAT IS COMPUTED
## A classical measurement-error reliability per cell,
##
##     lambda = (var(estimate) - mean(sampling variance)) / var(estimate),
##
## the fraction of between-locus variance in the estimate that is true
## between-locus variance rather than estimation noise. Under the usual
## attenuation argument, two independent estimates of the same underlying
## ordering correlate at sqrt(lambda_1 * lambda_2); that product is the ceiling
## the observed cross-cell correlation must be read against.
##
## WHAT THIS IS NOT
## This is an ANALYTIC ceiling, and it is an upper bound on the truth in two
## specific ways that must travel with the number:
##
##   * The per-locus uncertainty available is for the INPUT features
##     (bslmm_pve's posterior interval, he_h2's standard error), not for
##     pve_cis_joint_unbounded, which is what Stage 04 actually ranks. The
##     frozen joint model is a deterministic function of the features, so the
##     feature-level reliability transfers only to the extent the joint
##     prediction tracks that feature -- reported here per cell as
##     cor_with_score_basis so the reader can see how far the transfer is
##     being stretched.
##   * BSLMM's posterior SD is a shrunk Bayesian interval, not a frequentist
##     sampling SE, so lambda_bslmm is optimistic. he_se is frequentist but
##     HE's estimate is unbounded and its between-locus variance is inflated by
##     tails the score basis does not have, which pushes lambda_he optimistic
##     from the other direction. They are reported separately, never averaged.
##
## The confirmatory design is an empirical split-half: re-materialize each cell
## as two disjoint donor halves through 01b_estimation_cells, re-estimate, and
## correlate the two halves' within-half rank scores, Spearman-Brown adjusted
## back to full n. That measures the joint estimator end to end and needs no
## transfer argument. It is expensive -- at all_individuals.EA x dlpfc the
## halves are n=27 and n=28 -- and it is not run here. Until it is, the numbers
## below are a bound, not a measurement, and the manuscript should say so.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
})

opts <- parse_v2_args(require = c("region", "run_ids"))
region <- tolower(opts$region)
validate_cohort_region(region = region)

run_ids <- trimws(strsplit(opts$run_ids, ",", fixed = TRUE)[[1]])
run_ids <- run_ids[nzchar(run_ids)]
if (!length(run_ids) %in% c(1L, 2L)) {
    stop("Give one run ID (per-cell reliability) or two (reliability plus the ",
         "cross-cell ceiling): --run-ids a[,b]")
}

module_root <- file.path(V2_ROOT, "02_local_genetic_variance")

## Z for a central 95% interval: q975 - q025 spans 2 * 1.959964 posterior SDs.
Z95 <- 2 * stats::qnorm(0.975)

load_cell <- function(run_id) {
    run_dir <- file.path(module_root, "_m", "runs", run_id)
    if (!dir.exists(run_dir)) stop("Run directory not found: ", run_dir)
    man <- fread(file.path(run_dir, "manifest.tsv"))
    mval <- function(field) {
        v <- man$value[man$field == field]
        if (length(v) != 1L) stop("Run manifest lacks unique field: ", field)
        as.character(v[[1]])
    }
    if (!nzchar(mval("finished_at"))) {
        stop("Run is not sealed; a reliability ceiling from an open run is not ",
             "reproducible: ", run_id)
    }
    if (!identical(tolower(mval("region")), region)) {
        stop("Run is not for region ", region, ": ", run_id)
    }
    f <- list.files(file.path(run_dir, "results", "combined"),
                    pattern = "^local-genetic-control-.*-vmrs\\.tsv$",
                    full.names = TRUE)
    if (length(f) != 1L) {
        stop("Expected exactly one local-genetic-control table in ", run_id)
    }
    d <- fread(f[[1]])
    needed <- c("vmr_id", "local_genetic_control_eligible",
                "pve_cis_joint_unbounded", "bslmm_pve",
                "bslmm_pve_q025", "bslmm_pve_q975", "he_h2", "he_se")
    missing <- setdiff(needed, names(d))
    if (length(missing)) {
        stop("Score table for ", run_id, " lacks: ",
             paste(missing, collapse = ", "))
    }
    eligible <- tolower(trimws(as.character(d$local_genetic_control_eligible))) %in%
        c("true", "t", "1")
    d <- d[eligible]
    if (nrow(d) < 100L) {
        stop("Fewer than 100 eligible loci in ", run_id,
             "; reliability is not estimable")
    }
    list(run_id = run_id, cell = mval("cohort"), n = as.integer(mval("n_donors")),
         catalog_cohort = parse_cell(mval("cohort"))$catalog_cohort,
         vmr_set_id = mval("vmr_set_id"), d = d)
}

## lambda = 1 - mean(sampling variance) / var(estimate), clamped at 0. A
## negative value means the average per-locus uncertainty exceeds the total
## between-locus spread -- i.e. no recoverable ordering at all -- and is
## reported as 0 rather than as a negative correlation bound.
reliability <- function(estimate, se) {
    ok <- is.finite(estimate) & is.finite(se)
    estimate <- estimate[ok]
    se <- se[ok]
    if (length(estimate) < 100L) return(list(lambda = NA_real_, n = length(estimate)))
    total <- stats::var(estimate)
    noise <- mean(se^2)
    list(lambda = max(0, (total - noise) / total), n = length(estimate),
         var_total = total, var_noise = noise)
}

summarize_cell <- function(cell, keep_ids = NULL) {
    d <- cell$d
    if (!is.null(keep_ids)) d <- d[vmr_id %in% keep_ids]
    bslmm_sd <- (d$bslmm_pve_q975 - d$bslmm_pve_q025) / Z95
    lb <- reliability(d$bslmm_pve, bslmm_sd)
    lh <- reliability(d$he_h2, d$he_se)
    data.table(
        run_id = cell$run_id,
        cell = cell$cell,
        region = region,
        n_donors = cell$n,
        n_loci = nrow(d),
        lambda_bslmm = lb$lambda,
        lambda_he = lh$lambda,
        var_bslmm = lb$var_total,
        mean_sampling_var_bslmm = lb$var_noise,
        var_he = lh$var_total,
        mean_sampling_var_he = lh$var_noise,
        ## How far the feature-level reliability has to be carried to stand in
        ## for the reliability of the quantity Stage 04 actually ranks.
        cor_bslmm_with_score_basis = stats::cor(
            d$bslmm_pve, d$pve_cis_joint_unbounded,
            method = "spearman", use = "complete.obs"),
        cor_he_with_score_basis = stats::cor(
            d$he_h2, d$pve_cis_joint_unbounded,
            method = "spearman", use = "complete.obs"),
        ceiling_basis = "analytic_feature_level_upper_bound"
    )
}

cells <- lapply(run_ids, load_cell)

if (length(cells) == 2L) {
    if (identical(cells[[1]]$cell, cells[[2]]$cell)) {
        stop("Both runs are the same cell (", cells[[1]]$cell,
             "); a cross-cell ceiling needs two donor groups")
    }
    ## The donor-group contrast is two estimation groups over ONE discovery
    ## catalog (AGENTS.md 7.7). AA versus all_individuals is a set against its
    ## own superset, shares no vmr_set_id, and is not a donor-group contrast;
    ## refuse it here rather than emitting a ceiling for a quantity nobody may
    ## report.
    if (!identical(cells[[1]]$catalog_cohort, cells[[2]]$catalog_cohort)) {
        stop("Runs come from different discovery catalogs (",
             cells[[1]]$catalog_cohort, " vs ", cells[[2]]$catalog_cohort,
             "); a donor-group contrast shares one catalog (AGENTS.md 7.7)")
    }
    if (!identical(cells[[1]]$vmr_set_id, cells[[2]]$vmr_set_id)) {
        stop("Runs carry different vmr_set_id (", cells[[1]]$vmr_set_id, " vs ",
             cells[[2]]$vmr_set_id, "); the two cells must share one locus set")
    }
    shared <- intersect(cells[[1]]$d$vmr_id, cells[[2]]$d$vmr_id)
    if (length(shared) < 100L) {
        stop("Fewer than 100 loci eligible in both cells; the concordance ",
             "result has no usable support")
    }
    ## Reliability is computed on the SHARED eligible set, because that is the
    ## set the cross-cell correlation is computed on. A lambda from each cell's
    ## full eligible set would be the ceiling for a different quantity.
    out <- rbindlist(lapply(cells, summarize_cell, keep_ids = shared))
} else {
    shared <- NULL
    out <- summarize_cell(cells[[1]])
}

dest <- file.path(module_root, "_m", "combined",
                  sprintf("reliability-ceiling-%s.tsv", region))
dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

if (length(cells) == 2L) {
    d1 <- cells[[1]]$d[vmr_id %in% shared][order(vmr_id)]
    d2 <- cells[[2]]$d[vmr_id %in% shared][order(vmr_id)]
    stopifnot(identical(d1$vmr_id, d2$vmr_id))
    ## The observed quantity. Ranking within cell first is not decorative: the
    ## two cells have different eligible denominators, and AGENTS.md 7.6 forbids
    ## comparing raw scores across cells. Spearman on the raw basis and Spearman
    ## on the within-shared-set ranks are the same number; computing it this way
    ## makes the prohibition visible in the code.
    observed <- stats::cor(rank(d1$pve_cis_joint_unbounded),
                           rank(d2$pve_cis_joint_unbounded))
    ceil_b <- sqrt(out$lambda_bslmm[1] * out$lambda_bslmm[2])
    ceil_h <- sqrt(out$lambda_he[1] * out$lambda_he[2])
    out[, `:=`(
        n_shared_loci = length(shared),
        observed_cross_cell_spearman = observed,
        ceiling_bslmm = ceil_b,
        ceiling_he = ceil_h,
        ## Observed correlation as a fraction of what perfect agreement could
        ## have produced given this much estimation noise. This is the number
        ## the concordance claim rests on, NOT the raw Spearman.
        concordance_fraction_of_ceiling_bslmm = observed / ceil_b,
        concordance_fraction_of_ceiling_he = observed / ceil_h
    )]
}

write_tsv_plain <- function(x, path) {
    tmp <- tempfile(pattern = paste0(".", basename(path), "."),
                    tmpdir = dirname(path))
    write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    if (!file.rename(tmp, path)) stop("Could not atomically write: ", path)
    invisible(path)
}
write_tsv_plain(out, dest)

cat("\nReliability ceiling --", region, "\n")
for (i in seq_len(nrow(out))) {
    cat(sprintf("  %-26s n=%-4s loci=%-6s lambda_bslmm=%.3f lambda_he=%.3f",
                out$cell[i], out$n_donors[i], out$n_loci[i],
                out$lambda_bslmm[i], out$lambda_he[i]),
        sprintf(" (rho to score basis: bslmm %.3f, he %.3f)\n",
                out$cor_bslmm_with_score_basis[i],
                out$cor_he_with_score_basis[i]))
}
if (length(cells) == 2L) {
    cat(sprintf("\n  shared eligible loci        : %d\n", out$n_shared_loci[1]))
    cat(sprintf("  observed cross-cell Spearman: %.4f\n",
                out$observed_cross_cell_spearman[1]))
    cat(sprintf("  ceiling (bslmm)             : %.4f  -> %.1f%% of ceiling\n",
                out$ceiling_bslmm[1],
                100 * out$concordance_fraction_of_ceiling_bslmm[1]))
    cat(sprintf("  ceiling (he)                : %.4f  -> %.1f%% of ceiling\n",
                out$ceiling_he[1],
                100 * out$concordance_fraction_of_ceiling_he[1]))
    cat("\n  These ceilings are analytic upper bounds computed from the INPUT\n",
        " features' uncertainty, not from a re-estimation of the joint model.\n",
        " The empirical split-half confirmation has not been run.\n", sep = "")
}
cat("\nWrote", dest, "\n")
