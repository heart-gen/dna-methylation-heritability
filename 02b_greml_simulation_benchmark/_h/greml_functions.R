#### 02b_greml_simulation_benchmark -- shared functions ####
##
## Everything this module does to GCTA lives here, so arm 1 (the faithful v1
## port) and arm 2 (observed genotypes) call one implementation of the GRM
## writer, the REML call, the .hsq parser and the failure classifier.
##
## Nothing here reads or writes 02_local_genetic_variance/. The two observed-
## locus helpers this module needs -- load_observed_locus() and
## simulate_phenotype_on_observed_genotype() -- are in 00_shared/locus_io.R.

suppressPackageStartupMessages({
    library(data.table)
})

## ------------------------------------------------------------------ config
## Read config/greml_benchmark.yml from the run's code snapshot when one exists
## (see 09b_aging_application/_h/run_config.R for why), else the live tree.
load_greml_config <- function(run_dir = NULL) {
    if (!is.null(run_dir)) {
        f <- file.path(run_dir, "code", "config", "greml_benchmark.yml")
        if (file.exists(f)) {
            cfg <- yaml::read_yaml(f)
            attr(cfg, "config_file") <- f
            attr(cfg, "config_sha256") <- file_sha256(f)
            return(cfg)
        }
    }
    load_config("greml_benchmark")
}

## Refuse any config that would let this module say something it cannot.
assert_greml_interpretation <- function(cfg) {
    flags <- cfg$interpretation
    for (nm in c("absolute_pve_interpretation_allowed_for_observed_loci",
                 "observed_locus_estimation_allowed",
                 "heritability_class_allowed",
                 "reopens_module_02")) {
        if (!identical(flags[[nm]], FALSE)) {
            stop("config/greml_benchmark.yml interpretation.", nm,
                 " must be false (AGENTS.md 2.3, 3, 7.2).")
        }
    }
    if (!identical(flags$simulated_phenotypes_only, TRUE)) {
        stop("interpretation.simulated_phenotypes_only must be true: this ",
             "module never estimates an observed methylation phenotype.")
    }
    invisible(TRUE)
}

## Read a run manifest into a named character vector.
read_run_manifest <- function(run_dir) {
    m <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
    stats::setNames(m$value, m$field)
}

run_dir_for <- function(run_id) {
    file.path(repo_root(), "02b_greml_simulation_benchmark", "_m", "runs", run_id)
}

## The GCTA calls, the .hsq parser, the classifier and the GRM writer live in
## 00_shared/gcta.R (shared with 02c_cis_greml_sensitivity since 2026-10-03).

## ------------------------------------------------- arm 1: LD stratification
## Ported from simulation-analysis/gcta/_h/01.stratify_LD.R unchanged in
## substance: the breakpoints are summary()'s 1st quartile, median and 3rd
## quartile of the per-SNP LD score, and the four groups are
## (-Inf, Q1], (Q1, median], (median, Q3], (Q3, Inf).
stratify_ld_quartiles <- function(ldscore) {
    q <- summary(ldscore)
    g <- integer(length(ldscore))
    g[ldscore <= q[[2]]] <- 1L
    g[ldscore > q[[2]] & ldscore <= q[[3]]] <- 2L
    g[ldscore > q[[3]] & ldscore <= q[[5]]] <- 3L
    g[ldscore > q[[5]]] <- 4L
    list(group = g, breaks = c(q1 = q[[2]], median = q[[3]], q3 = q[[5]]))
}

## ------------------------------------------------------------------ metrics
## Continuous recovery metrics against the realized simulated h2. No class, no
## threshold: AGENTS.md 2.3 and 3 retire both.
recovery_metrics <- function(est, se, truth, constrained) {
    ok <- is.finite(est) & is.finite(truth)
    est <- est[ok]; se <- se[ok]; truth <- truth[ok]
    n <- length(est)
    if (n == 0L) {
        return(data.table(n_estimates = 0L, mean_truth = NA_real_,
                          mean_estimate = NA_real_, bias = NA_real_,
                          rmse = NA_real_, coverage95 = NA_real_,
                          boundary_rate = NA_real_,
                          outside_unit_interval_rate = NA_real_))
    }
    covered <- is.finite(se) & abs(est - truth) <= stats::qnorm(0.975) * se
    data.table(
        n_estimates = n,
        mean_truth = mean(truth),
        mean_estimate = mean(est),
        bias = mean(est - truth),
        rmse = sqrt(mean((est - truth)^2)),
        coverage95 = mean(covered),
        ## Constrained REML piles estimates onto the [0, 1] bounds; the
        ## unconstrained estimator instead leaves the interval. Report the one
        ## that describes each mode, and the other as NA.
        boundary_rate = if (constrained)
            mean(est <= 1e-5 | est >= 1 - 1e-5) else NA_real_,
        outside_unit_interval_rate = if (!constrained)
            mean(est < 0 | est > 1) else NA_real_
    )
}

## Spearman between truth and estimate with a cluster bootstrap CI. Clusters
## are loci in arm 2 (scenarios on one locus share its genotype) and
## phenotypes in arm 1.
spearman_cluster_ci <- function(est, truth, cluster, n_boot, seed) {
    ok <- is.finite(est) & is.finite(truth)
    est <- est[ok]; truth <- truth[ok]; cluster <- cluster[ok]
    if (length(est) < 3L) {
        return(data.table(spearman = NA_real_, ci_low = NA_real_,
                          ci_high = NA_real_, n_clusters = length(unique(cluster))))
    }
    rho <- suppressWarnings(stats::cor(est, truth, method = "spearman"))
    cl <- unique(cluster)
    idx <- split(seq_along(cluster), cluster)
    set.seed(seed)
    boots <- vapply(seq_len(n_boot), function(b) {
        pick <- unlist(idx[sample(as.character(cl), length(cl), replace = TRUE)],
                       use.names = FALSE)
        suppressWarnings(stats::cor(est[pick], truth[pick], method = "spearman"))
    }, numeric(1))
    ci <- stats::quantile(boots, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    data.table(spearman = rho, ci_low = ci[[1]], ci_high = ci[[2]],
               n_clusters = length(cl))
}

## The one definition of an arm 2 unit ID, used by the stage that writes a
## unit's status and by the stage that reconciles it. Two inline sprintf()s
## disagreed on the first smoke run: format() over the whole h2 vector pads to
## a common width ("0.00"), over one value it does not ("0"), so 144 of 192
## units looked unexpected. as.character() formats element-wise.
arm2_unit_id <- function(vmr_id, h2, architecture, replicate, reml_mode) {
    sprintf("%s|h2=%s|%s|rep%d|%s", vmr_id, as.character(h2), architecture,
            as.integer(replicate), reml_mode)
}

## Write the per-task status rows every stage emits, so reconciliation counts
## units rather than trusting a SLURM exit code.
write_status <- function(rows, path) {
    rows <- as.data.table(rows)
    stopifnot(all(c("unit", "status", "reason") %in% names(rows)))
    write_atomic(rows, path)
}
