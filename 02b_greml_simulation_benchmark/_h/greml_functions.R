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

## ------------------------------------------------------------------- GCTA
gcta_version <- function(bin) {
    out <- suppressWarnings(system2(bin, character(0), stdout = TRUE,
                                    stderr = TRUE))
    hit <- regmatches(out, regexpr("v[0-9]+\\.[0-9]+(\\.[0-9]+)?", out))
    if (!length(hit)) stop("Could not read a version from ", bin)
    sub("^v", "", hit[[1]])
}

## Stop unless the binary is the pinned one. A REML result depends on the
## version (defaults and convergence criteria have changed across releases).
assert_gcta_pin <- function(cfg) {
    bin <- cfg$gcta$binary
    if (!file.exists(bin)) stop("GCTA binary not found: ", bin)
    v <- gcta_version(bin)
    if (!identical(v, as.character(cfg$gcta$version))) {
        stop("GCTA version ", v, " does not match the pin ",
             cfg$gcta$version, " in config/greml_benchmark.yml")
    }
    v
}

## Run gcta64 and classify the outcome.
##
## status is "completed" when an .hsq was written; "qc_failed" when GCTA stopped
## with one of the configured REML estimation messages (singular information
## matrix, non-convergence) -- an outcome of the estimator on that phenotype,
## counted and reported, not a crash; and "failed" for anything else, which
## reconcile() treats as an unexplained computational failure.
run_gcta <- function(bin, args, out_prefix, threads = 1L) {
    log <- paste0(out_prefix, ".gcta.log")
    st <- system2(bin, c(args, "--out", out_prefix, "--thread-num",
                         as.character(threads)),
                  stdout = log, stderr = log)
    invisible(list(exit = st, log = log))
}

classify_reml <- function(exit, hsq, log, failure_patterns) {
    if (identical(as.integer(exit), 0L) && file.exists(hsq)) {
        return(list(status = "completed", reason = NA_character_))
    }
    txt <- if (file.exists(log)) paste(readLines(log, warn = FALSE),
                                       collapse = "\n") else ""
    for (p in failure_patterns) {
        if (grepl(p, txt, fixed = TRUE)) {
            return(list(status = "qc_failed",
                        reason = paste0("reml_estimation_failure: ", p)))
        }
    }
    tail_txt <- utils::tail(strsplit(txt, "\n")[[1]], 3)
    list(status = "failed",
         reason = paste0("gcta_exit_", exit, ": ",
                         paste(tail_txt, collapse = " | ")))
}

## Parse a GCTA .hsq into one row. Handles single-GRM ("V(G)/Vp") and
## multi-GRM ("Sum of V(G)/Vp") output; the multi-GRM total is what v1's
## 02.summary.py reported as `Sum of V(G)_Vp_Variance`.
parse_hsq <- function(path) {
    lines <- readLines(path, warn = FALSE)
    parts <- strsplit(lines[-1], "\t", fixed = TRUE)
    src <- vapply(parts, `[`, character(1), 1)
    val <- suppressWarnings(as.numeric(vapply(parts, `[`, character(1), 2)))
    se <- suppressWarnings(as.numeric(vapply(parts, function(x)
        if (length(x) >= 3) x[[3]] else NA_character_, character(1))))
    get <- function(s, v = val) { i <- match(s, src); if (is.na(i)) NA_real_ else v[[i]] }
    total <- if ("Sum of V(G)/Vp" %in% src) "Sum of V(G)/Vp" else "V(G)/Vp"
    data.table(
        h2_hat = get(total), h2_se = get(total, se),
        vg_hat = get("V(G)"), ve_hat = get("V(e)"), vp_hat = get("Vp"),
        logL = get("logL"), logL0 = get("logL0"), lrt = get("LRT"),
        pval = get("Pval"), n_used = get("n")
    )
}

## ------------------------------------------------------------------ the GRM
## A GCTA-format GRM written directly from a dosage matrix. Matches
## `gcta64 --make-grm` on the same SNPs and donors to float32 precision
## (max |diff| 1.0e-7 on a 118-donor, 3,428-SNP DLPFC locus with 1,930 missing
## calls; tests/test_greml_functions.R repeats that check). Missing calls
## contribute nothing and the denominator is the per-pair count of SNPs
## observed in both donors, which is GCTA's own rule. Writing it here rather
## than re-exporting a PLINK file keeps the genotype exactly the matrix that
## load_observed_locus() returned -- same window, same QC, same donor order.
grm_from_dosage <- function(G) {
    p <- colMeans(G, na.rm = TRUE) / 2
    keep <- is.finite(p) & p > 0 & p < 1
    G <- G[, keep, drop = FALSE]; p <- p[keep]
    Z <- sweep(G, 2, 2 * p, "-")
    Z <- sweep(Z, 2, sqrt(2 * p * (1 - p)), "/")
    M <- !is.na(Z)
    Z[!M] <- 0
    N <- tcrossprod(M * 1)
    list(A = tcrossprod(Z) / N, N = N, n_snps = ncol(G))
}

## GCTA stores the lower triangle row by row, which is the upper triangle in
## R's column-major order.
write_grm <- function(grm, ids, prefix) {
    A <- grm$A; N <- grm$N
    writeBin(as.numeric(A[upper.tri(A, diag = TRUE)]),
             paste0(prefix, ".grm.bin"), size = 4)
    writeBin(as.numeric(N[upper.tri(N, diag = TRUE)]),
             paste0(prefix, ".grm.N.bin"), size = 4)
    fwrite(as.data.table(ids)[, 1:2], paste0(prefix, ".grm.id"),
           sep = "\t", col.names = FALSE)
    invisible(prefix)
}

read_grm <- function(prefix) {
    ids <- fread(paste0(prefix, ".grm.id"), header = FALSE)
    n <- nrow(ids)
    v <- readBin(paste0(prefix, ".grm.bin"), "numeric", size = 4,
                 n = n * (n + 1) / 2)
    A <- matrix(0, n, n)
    A[upper.tri(A, diag = TRUE)] <- v
    A <- A + t(A) - diag(diag(A))
    list(A = A, ids = ids)
}

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
