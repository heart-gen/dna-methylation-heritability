#!/usr/bin/env Rscript
#### 07 Stage 08 -- power of each modality's coupling tests ####
##
## Usage:
##   Rscript _h/08_power_analysis.R --run-id tsc-AA-<region>-YYYYMMDD
##   Rscript _h/08_power_analysis.R --run-id <id> --out-dir <dir>   # read-only
##
## Why this stage exists
## --------------------
## `config/transcription_splicing.yml:gates:min_vmrs_tested` is a locked power
## floor, and 03_test_coupling.R keeps a modality below it out of the
## coupling-test FDR family. A floor stated as a VMR count is easy to read as
## arbitrary, so this stage states the same restriction as the quantity a reader
## actually cares about: the smallest association each modality could have
## detected, given the design it actually had.
##
## The binary predictor is the case that forces the issue. `any_meqtl_support`
## crossed with the coupled/not-coupled outcome is a 2x2, and on 2026-09-25 its
## off-cell was EMPTY in all three AA regions -- 8 to 15 coupled VMRs, every one
## of them with meQTL support. A logistic fit on that table has no maximum
## likelihood estimate: the coefficient runs to the boundary, and the reported
## estimate near 17.6 with an SE near 0.4 and p underflowing to 0 is an
## optimizer stopping, not an effect. Reporting the minimum detectable odds ratio
## makes clear that this is a property of the ABC link set's SIZE and not of the
## data it happened to contain.
##
## Method
## ------
## For the binary predictor, planning is done on the 2x2 with the OBSERVED
## margins and cells expected under independence, which is the standard
## fixed-margin calculation and does not depend on the observed association:
##
##   SE(log OR) = sqrt(sum_ij 1 / E_ij),   E_ij = row_i * col_j / N
##   minimum detectable |log OR| at power 1-beta, two-sided alpha
##            = (z_{1-alpha/2} + z_{1-beta}) * SE(log OR)
##
## Using expected rather than observed cells is deliberate: an observed zero cell
## makes SE infinite, so an observed-cell calculation would report "no detectable
## effect" for the very design whose power we are trying to quantify, and would
## also make the answer depend on the result. The margins are design facts; the
## cross-classification is the result.
##
## For the continuous predictors the same z-based expression is used with
## SE(log OR) approximated from the outcome prevalence, which is the usual
## logistic planning approximation:
##
##   SE(beta) ~ 1 / sqrt(N * p * (1 - p) * var(x))
##
## Both are planning approximations and are labelled as such in the output. They
## are not a substitute for the fitted SEs, which the coupling-tests table
## already carries; they answer a different question -- what this design could
## have found -- and only that.
##
## Every row carries `separated`, so a reader can see that the exclusion and the
## separation are two consequences of one cause rather than two findings.

V2_TSC_H <- Sys.getenv("V2_RUN_CODE", file.path(Sys.getenv("V2_REPO_ROOT", "."),
                                                "07_transcription_splicing_coupling", "_h"))
source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "07_transcription_splicing_coupling"
opts <- parse_v2_args(require = c("run_id"))
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

ts <- load_config("transcription_splicing")
ALPHA <- 0.05
POWER <- 0.80
Z <- stats::qnorm(1 - ALPHA / 2) + stats::qnorm(POWER)

mf_f <- file.path(run_dir, "results", "coupling-model-frame.tsv")
if (!file.exists(mf_f)) {
    stop("No coupling-model-frame.tsv under ", opts$run_id,
         "; run _h/03_test_coupling.R first.")
}
mf <- fread(mf_f)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]
    if (length(v) != 1L) NA_character_ else as.character(v[[1L]])
}

## ------------------------------------------------------------------ binary arm
binary_row <- function(d, mod) {
    tb <- table(factor(d$coupled, levels = c(0, 1)),
                factor(as.integer(d$any_meqtl_support > 0), levels = c(0, 1)))
    N <- sum(tb)
    rs <- rowSums(tb); cs <- colSums(tb)
    E <- outer(rs, cs) / N
    ## A zero EXPECTED cell means a margin is degenerate -- no coupled VMRs at
    ## all, or none with meQTL support -- and then nothing is estimable by any
    ## method. Distinguish that from a zero OBSERVED cell, which is separation.
    degenerate <- any(E == 0)
    se_log_or <- if (degenerate) NA_real_ else sqrt(sum(1 / E))
    data.table(
        modality = mod,
        predictor = "any_meqtl_support",
        predictor_type = "binary",
        n_vmrs = N,
        n_coupled = as.integer(rs[["1"]]),
        n_with_meqtl_support = as.integer(cs[["1"]]),
        prevalence_coupled = as.numeric(rs[["1"]]) / N,
        min_observed_cell = as.integer(min(tb)),
        min_expected_cell = round(min(E), 2),
        separated = min(tb) == 0,
        degenerate_margin = degenerate,
        se_log_or_expected = round(se_log_or, 4),
        min_detectable_or = round(exp(Z * se_log_or), 2),
        method = "2x2 fixed-margin normal approximation on log OR")
}

## -------------------------------------------------------------- continuous arm
continuous_row <- function(d, mod, pcol, pname) {
    x <- d[[pcol]]
    keep <- is.finite(x) & is.finite(d$coupled)
    x <- x[keep]; y <- d$coupled[keep]
    N <- length(y); p <- mean(y)
    vx <- stats::var(x)
    se_beta <- if (N == 0 || p <= 0 || p >= 1 || !is.finite(vx) || vx <= 0) {
        NA_real_
    } else {
        1 / sqrt(N * p * (1 - p) * vx)
    }
    data.table(
        modality = mod,
        predictor = pname,
        predictor_type = "continuous",
        n_vmrs = N,
        n_coupled = sum(y),
        n_with_meqtl_support = NA_integer_,
        prevalence_coupled = p,
        min_observed_cell = NA_integer_,
        min_expected_cell = NA_real_,
        separated = FALSE,
        degenerate_margin = !is.finite(se_beta),
        se_log_or_expected = round(se_beta, 4),
        min_detectable_or = round(exp(Z * se_beta), 2),
        method = "logistic planning approximation, SE ~ 1/sqrt(N p (1-p) var(x))")
}

rows <- list()
for (mod in unique(mf$modality)) {
    d <- mf[modality == mod]
    if (!nrow(d)) next
    if ("any_meqtl_support" %in% names(d)) {
        rows[[length(rows) + 1]] <- binary_row(d, mod)
    }
    for (spec in list(
            c("proportion_cpgs_with_sig_meqtl", "meqtl_proportion"),
            c("local_snp_contribution_score_z", "local_genetic_control"))) {
        if (spec[[1]] %in% names(d)) {
            rows[[length(rows) + 1]] <- continuous_row(d, mod, spec[[1]], spec[[2]])
        }
    }
}
if (!length(rows)) stop("No modality produced a power row")
out <- rbindlist(rows, fill = TRUE)

floor_vmrs <- ts$gates$min_vmrs_tested
out[, power_floor_min_vmrs := floor_vmrs]
out[, meets_power_floor := n_vmrs >= floor_vmrs]
out[, in_fdr_family := meets_power_floor]
out[, `:=`(alpha = ALPHA, target_power = POWER,
           cohort = mval("cohort"), region = mval("region"),
           run_id = opts$run_id, vmr_set_id = mval("vmr_set_id"))]
out[, planning_estimate_not_a_fitted_se := TRUE]
setorder(out, modality, predictor)

dest_dir <- if (!is.null(opts$out_dir)) opts$out_dir else
    file.path(run_dir, "results")
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
dest <- file.path(dest_dir, "coupling-power-analysis.tsv")
if (is.null(opts$out_dir)) {
    write_atomic(out, dest)
} else {
    fwrite(out, dest, sep = "\t")
}

print(out[, .(modality, predictor, n_vmrs, n_coupled, min_observed_cell,
              separated, min_detectable_or, meets_power_floor)])
message("[07] power analysis -> ", dest)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
