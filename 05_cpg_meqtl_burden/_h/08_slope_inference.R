#!/usr/bin/env Rscript
#### 05_cpg_meqtl_burden -- donor-robust inference for the DLPFC-hippocampus slope ####
##
## Usage:
##   Rscript _h/08_slope_inference.R --run-id cmb-AA-crossregion-YYYYMMDD
##
## THE QUESTION. Module 08 tier 2 reads the CpG meQTL-burden slope as differing
## between DLPFC (2.094) and hippocampus (2.367): z -3.63, q 0.0038. That z uses
## sqrt(se_a^2 + se_b^2) over 03_vmr_burden.R's VMR-level HC3 SEs. Two things
## are wrong with it, in opposite directions, so its calibration is unknown:
##   - HC3 treats ~9,000 VMRs as independent. Nearby VMRs are not, and every VMR
##     in a region is estimated in the same donors (anti-conservative).
##   - The independence form ignores that the two regions share 115 of 118
##     donors, which correlates the two slopes positively (conservative).
## AGENTS.md 7.5 requires donor-robust inference. This stage supplies it, and
## until it existed the writing rule was: report the slopes, claim no
## heterogeneity (TASKS.md A3).
##
## THE VARIANCE. The construction 00_shared/axis_inference.R uses and T28
## validated for 09b and 10:
##   SE^2 = paired donor-bootstrap variance + JOINT delete-one-chromosome
##          weighted block-jackknife variance
## Donor half: 07_bootstrap_map.py re-mapped every tested CpG on each paired
## draw and 08a_bootstrap_counts.py re-counted the burden, so each draw refits
## this model on a re-derived OUTCOME with the predictor and VMR covariates held
## fixed. VMR half: chromosome k removed from both regions at once, on the
## accepted tables.
##
## WHAT T28 DID NOT VALIDATE. T28 calibrated this sum for a debiased
## squared-effect outcome. This outcome is a count of significant CpGs, and a
## bootstrap draw with duplicated donors calls more CpGs significant than the
## original sample did (+18% at p < 1e-4 on the chr22 probe), so the bootstrap
## distribution is shifted. The bootstrap is used for its VARIANCE only, never
## for a percentile interval, and the shift is reported (bootstrap_mean_*), not
## hidden. Calibration of the sum for this outcome has not been simulated.
##
## Nothing here claims heterogeneity. Whether a difference is claimed is Module
## 08 tier 2's decision under its own strict conjunction; this run supplies the
## variance that decision should use.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(sandwich)
    library(lmtest)
})

MODULE <- "05_cpg_meqtl_burden"
PRED <- "local_snp_contribution_score_z"
REPRO_TOL <- 1e-6

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]; if (length(v) == 0) NA_character_ else v[1]
}
smoke <- identical(mval("smoke_run"), "TRUE")
contrast <- strsplit(mval("contrast"), ",")[[1]]
B <- as.integer(mval("bootstrap_n"))
si <- load_config("meqtl_parameters")$cross_region_slope_inference
z_ci <- stats::qnorm(1 - (1 - as.numeric(si$ci_level)) / 2)
res_dir <- file.path(run_dir, "results")

## ------------------------------------------- the accepted model, refitted
cell_run <- function(re) mval(paste0("upstream_cpg_meqtl_burden_", re))
cell_res <- function(re) file.path(repo_root(), MODULE, "_m", "runs", cell_run(re), "results")

FORM <- cbind(n_sig, n_tested_cpgs - n_sig) ~ local_snp_contribution_score_z +
    log(vmr_length) + cpg_density + mean_methylation

slope <- function(d) {
    fit <- tryCatch(suppressWarnings(stats::glm(FORM, data = d,
                                                family = stats::quasibinomial())),
                    error = function(e) NULL)
    if (is.null(fit)) NA_real_ else unname(stats::coef(fit)[PRED])
}

burden <- lapply(contrast, function(re) {
    d <- fread(file.path(cell_res(re), "vmr-meqtl-burden.tsv"))
    d[, n_sig := n_cpgs_with_sig_meqtl]
    d[]
})
names(burden) <- contrast

## Reproduction guard: the slope and HC3 SE refitted here must be the sealed
## ones, or the variance below belongs to a different model.
repro <- rbindlist(lapply(contrast, function(re) {
    d <- burden[[re]]
    fit <- stats::glm(FORM, data = d, family = stats::quasibinomial())
    ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
    s <- fread(file.path(cell_res(re), "burden-primary-model.tsv"))[term == PRED]
    data.table(region = re, run_id = cell_run(re), n_vmrs = nrow(d),
               estimate_refit = ct[PRED, 1], estimate_sealed = s$estimate,
               se_hc3_refit = ct[PRED, 2], se_hc3_sealed = s$se)
}))
if (repro[, any(abs(estimate_refit - estimate_sealed) > REPRO_TOL |
                abs(se_hc3_refit - se_hc3_sealed) > REPRO_TOL)]) {
    print(repro)
    stop("The refitted burden model does not reproduce the accepted cells")
}
obs <- stats::setNames(repro$estimate_refit, contrast)
d_obs <- obs[[1]] - obs[[2]]

## ------------------------------------------------- VMR half: chromosome jackknife
chroms <- sort(unique(unlist(lapply(burden, `[[`, "chrom"))))
jk <- rbindlist(lapply(chroms, function(k) {
    est <- vapply(contrast, function(re) slope(burden[[re]][chrom != k]), numeric(1))
    data.table(chrom = k,
               n_a = burden[[1]][chrom == k, .N], n_b = burden[[2]][chrom == k, .N],
               estimate_a = est[[1]], estimate_b = est[[2]],
               difference = est[[1]] - est[[2]])
}))
## Busing et al. weighted delete-a-block jackknife, as in 09b stage 05.
jk_var <- function(full, theta, n_block, n_total) {
    h <- n_total / n_block
    pseudo <- h * full - (h - 1) * theta
    sum((pseudo - mean(pseudo))^2 / (h - 1)) / length(theta)
}
var_jk_a <- jk_var(obs[[1]], jk$estimate_a, jk$n_a, sum(jk$n_a))
var_jk_b <- jk_var(obs[[2]], jk$estimate_b, jk$n_b, sum(jk$n_b))
var_jk_d <- jk_var(d_obs, jk$difference, jk$n_a + jk$n_b, sum(jk$n_a + jk$n_b))

## ------------------------------------------------ donor half: paired bootstrap
counts <- fread(file.path(res_dir, "bootstrap-burden-counts.tsv.gz"))
boot <- rbindlist(lapply(contrast, function(re) {
    cr <- counts[region == re]
    ## Draw 0 is the accepted cell pushed through the counting code; 08a has
    ## already checked it equals the sealed counts. Check the join here too.
    d0 <- merge(burden[[re]][, .(vmr_id, n_sealed = n_sig)], cr[draw == 0L],
                by = "vmr_id")
    if (nrow(d0) != nrow(burden[[re]]) || d0[, any(n_sealed != n_sig)]) {
        stop(re, ": draw 0 does not match the sealed burden table")
    }
    base <- burden[[re]][, !"n_sig"]
    rbindlist(lapply(seq_len(B), function(b) {
        d <- merge(base, cr[draw == b, .(vmr_id, n_sig)], by = "vmr_id")
        if (nrow(d) != nrow(base)) stop(re, " draw ", b, ": VMR set changed")
        data.table(region = re, draw = b, estimate = slope(d))
    }))
}))
bw <- dcast(boot, draw ~ region, value.var = "estimate")
setcolorder(bw, c("draw", contrast))
ok <- stats::complete.cases(bw)
ba <- bw[[contrast[1]]][ok]; bb <- bw[[contrast[2]]][ok]
var_boot_a <- stats::var(ba); var_boot_b <- stats::var(bb)
var_boot_d <- stats::var(ba - bb)

## ------------------------------------------------------------------ results
per_region <- data.table(
    region = contrast, run_id = vapply(contrast, cell_run, character(1)),
    n_vmrs = repro$n_vmrs, estimate = unname(obs),
    se_donor_robust = sqrt(c(var_boot_a + var_jk_a, var_boot_b + var_jk_b)),
    se_bootstrap = sqrt(c(var_boot_a, var_boot_b)),
    se_jackknife = sqrt(c(var_jk_a, var_jk_b)),
    se_hc3_sealed = repro$se_hc3_sealed,
    bootstrap_mean = c(mean(ba), mean(bb)))
per_region[, `:=`(ci_lower = estimate - z_ci * se_donor_robust,
                  ci_upper = estimate + z_ci * se_donor_robust,
                  se_ratio_to_hc3 = se_donor_robust / se_hc3_sealed)]

se_d <- sqrt(var_boot_d + var_jk_d)
difference <- data.table(
    region_a = contrast[1], region_b = contrast[2],
    estimate_a = obs[[1]], estimate_b = obs[[2]], difference = d_obs,
    se = se_d, se_bootstrap = sqrt(var_boot_d), se_jackknife = sqrt(var_jk_d),
    z = d_obs / se_d, p = 2 * stats::pnorm(-abs(d_obs / se_d)),
    ci_lower = d_obs - z_ci * se_d, ci_upper = d_obs + z_ci * se_d,
    ci_excludes_zero = (d_obs - z_ci * se_d > 0) | (d_obs + z_ci * se_d < 0),
    se_module08_independent_hc3 = sqrt(sum(repro$se_hc3_sealed^2)),
    z_module08_independent_hc3 = d_obs / sqrt(sum(repro$se_hc3_sealed^2)),
    bootstrap_cor_between_regions = stats::cor(ba, bb),
    bootstrap_mean_difference = mean(ba - bb),
    n_bootstrap = B, n_bootstrap_failed = sum(!ok),
    bootstrap_nperm = as.integer(mval("bootstrap_nperm")),
    n_chromosome_blocks = length(chroms),
    n_union_donors = as.integer(mval("n_union_donors")),
    n_shared_donors = as.integer(mval("n_shared_donors")),
    inference = "paired_donor_bootstrap_var_plus_joint_chromosome_jackknife_var",
    heterogeneity_claim = "not_made_here_module08_tier2_decides")

write_atomic(repro, file.path(res_dir, "slope-reproduction-check.tsv"))
write_atomic(per_region, file.path(res_dir, "slope-per-region.tsv"))
write_atomic(difference, file.path(res_dir, "slope-difference.tsv"))
write_atomic(jk, file.path(res_dir, "slope-jackknife.tsv"))
write_atomic(bw, file.path(res_dir, "slope-bootstrap-draws.tsv"))
writeLines(c(
    "Interpretation constraints carried by this run:",
    "  - Variance = paired donor-bootstrap + joint chromosome-jackknife, the",
    "    00_shared/axis_inference.R construction. T28 validated it for a",
    "    debiased squared-effect outcome, NOT for this significant-CpG count;",
    "    its calibration here is unsimulated.",
    "  - Bootstrap draws re-map meQTLs with duplicated donors and call more CpGs",
    "    significant than the original sample, so the bootstrap distribution is",
    "    shifted (bootstrap_mean_*). It is used for its variance only.",
    sprintf("  - Bootstrap mapping used %s permutations (the accepted cells used",
            mval("bootstrap_nperm")),
    "    tensorqtl's default 10,000); see config/meqtl_parameters.yml.",
    "  - Predictor, VMR covariates and the M3a covariate values are held fixed:",
    "    the variance is conditional on the locked design.",
    "  - DLPFC and hippocampus share sequencing batches 1-2, so this contrast is",
    "    the identified one (AGENTS.md 7.7 tier 2). No heterogeneity claim is",
    "    made here; Module 08 tier 2 decides."
), file.path(res_dir, "interpretation-constraints.txt"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "session-info.txt"))

print(per_region, digits = 4)
print(difference[, .(difference, se, se_bootstrap, se_jackknife, z, p,
                     ci_lower, ci_upper, z_module08_independent_hc3,
                     bootstrap_cor_between_regions)], digits = 4)

if (!smoke) {
    append_manifest(list(dir = run_dir), list(
        decision = if (difference$ci_excludes_zero) "SLOPE_DIFFERENCE_CI_EXCLUDES_ZERO"
                   else "SLOPE_DIFFERENCE_CI_INCLUDES_ZERO",
        sealed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
    close_run(list(dir = run_dir))
}
