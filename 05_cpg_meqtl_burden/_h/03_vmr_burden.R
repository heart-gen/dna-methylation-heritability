#!/usr/bin/env Rscript
#### 05_cpg_meqtl_burden -- VMR-level burden against relative local control ####
##
## Usage:
##   Rscript _h/03_vmr_burden.R --run-id cmb-AA-caudate-20260817
##
## The orthogonal-evidence test: does a VMR with a higher relative local SNP
## contribution contain a greater FRACTION of CpGs with conventional cis-meQTL
## support?
##
## Three design commitments, all from AGENTS.md 7.5:
##
##  1. The predictor is the CONTINUOUS local SNP contribution score. The
##     top-versus-bottom-quartile
##     (high vs low) comparison is secondary evidence and is computed second, in
##     its own table, so it cannot be mistaken for the primary result.
##  2. The outcome is a count of significant CpGs out of TESTED CpGs -- a
##     binomial proportion with a known denominator, so it is modelled as such
##     with overdispersion allowed, not as a continuous fraction.
##  3. This is convergent evidence from the SAME donors, not independent
##     replication, and the output says so.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(sandwich)
    library(lmtest)
})

MODULE <- "05_cpg_meqtl_burden"

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
res_dir <- file.path(run_dir, "results")

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]; if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mval("cohort"); region <- mval("region")
meqtl <- load_config("meqtl_parameters")

## ------------------------------------------------------------------- inputs
cpg_res <- fread(file.path(res_dir, "cpg-meqtl-results.tsv"))
denom   <- fread(file.path(res_dir, "vmr-cpg-denominators.tsv"))
member  <- fread(file.path(res_dir, "tested-cpg-membership.tsv"))

lcg <- load_local_genetic_control(
    mval("upstream_local_genetic_variance_run_id"),
    region = region, cohort = cohort, eligible_only = TRUE
)

## ------------------------------------------------------------ inflation gate
## AGENTS.md 7.5: resolve genomic inflation BEFORE the figure freeze. A burden
## computed on an inflated scan is inflated everywhere, and uniformly enough
## that it still produces a clean-looking gradient.
## Lambda comes from the NOMINAL pass, computed by 04_qc_plots.py, which runs
## before this stage. It cannot be recomputed from `cpg_res`: that table holds
## the single top variant per CpG, selected for extremeness, so its p-value
## distribution is not the null and its median chi-square is not an inflation
## estimate.
infl_f <- file.path(res_dir, "qc", "genomic-inflation.tsv")
if (!file.exists(infl_f)) {
    stop("Missing ", infl_f, ". Run 04_qc_plots.py before 03_vmr_burden.R; ",
         "the inflation gate is computed there from the nominal cis pairs.")
}
lambda <- fread(infl_f)$lambda[1]
lambda_max <- config_get(meqtl, "genomic_inflation.max")
message("[05] genomic inflation lambda = ", round(lambda, 3),
        " (distal cis pairs; allowed max ", lambda_max, ")")
if (!is.finite(lambda) || lambda > lambda_max) {
    stop("Genomic inflation lambda = ", round(lambda, 3),
         " exceeds the configured maximum of ", lambda_max,
         ". Resolve the inflation (covariates, relatedness, population ",
         "structure), or revise config/meqtl_parameters.yml:genomic_inflation ",
         "-- see the caveat recorded there about lambda over cis pairs.")
}

## --------------------------------------------------------------- aggregation
sig <- cpg_res[qvalue <= meqtl$mapping$fdr_threshold]
burden <- member[, .(vmr_id = unique(vmr_id)), by = cpg_id][
    , .(n_sig = sum(cpg_id %in% sig$cpg_id)), by = vmr_id]

burden <- merge(denom, burden, by = "vmr_id", all.x = TRUE)
burden[is.na(n_sig), n_sig := 0L]

## Denominator audit: a VMR with zero tested CpGs has no defined burden and is
## excluded here rather than contributing a 0/0 that R would silently make NaN.
untestable <- burden[n_tested_cpgs == 0]
burden <- burden[n_tested_cpgs > 0]
write_atomic(untestable, file.path(res_dir, "vmrs-without-tested-cpgs.tsv"))

burden <- merge(burden,
                lcg[, .(vmr_id, local_snp_contribution_score,
                        local_snp_contribution_score_z,
                        local_snp_contribution_quartile,
                        chrom, start, end, n_cpgs, mean_methylation)],
                by = "vmr_id")
burden[, `:=`(
    prop_sig = n_sig / n_tested_cpgs,
    vmr_length = end - start + 1L
)]
burden[, cpg_density := n_tested_cpgs / vmr_length]

## ------------------------------------------------------------ primary model
## Quasibinomial on the count/denominator pair: the denominator varies by an
## order of magnitude across VMRs, and treating prop_sig as a plain continuous
## outcome throws that precision information away.
fit <- stats::glm(cbind(n_sig, n_tested_cpgs - n_sig) ~
                      local_snp_contribution_score_z + log(vmr_length) +
                      cpg_density + mean_methylation,
                  data = burden, family = stats::quasibinomial())
ct <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))

primary <- data.table(
    model = "quasibinomial_continuous_local_snp_contribution_score",
    term = rownames(ct), estimate = ct[, 1], se = ct[, 2],
    z = ct[, 3], p = ct[, 4],
    n_vmrs = nrow(burden),
    dispersion = summary(fit)$dispersion,
    genomic_inflation_lambda = lambda,
    region = region, population = cohort, vmr_set_id = mval("vmr_set_id")
)
write_atomic(primary, file.path(res_dir, "burden-primary-model.tsv"))

## ------------------------------------------- secondary: matched extreme groups
## Secondary evidence only. Matched on the same covariates the primary model
## adjusts for, so the comparison is not just a restatement of VMR size.
ext <- burden[local_snp_contribution_quartile %in%
              c("bottom_quartile", "top_quartile")]
ext[, group := local_snp_contribution_quartile]
secondary <- if (nrow(ext) >= 40 && uniqueN(ext$group) == 2) {
    tt <- stats::wilcox.test(prop_sig ~ group, data = ext)
    data.table(model = "top_bottom_quartiles_SECONDARY",
               n_high = ext[group == "top_quartile", .N],
               n_low = ext[group == "bottom_quartile", .N],
               median_prop_sig_high = stats::median(
                   ext[group == "top_quartile"]$prop_sig),
               median_prop_sig_low = stats::median(
                   ext[group == "bottom_quartile"]$prop_sig),
               p = tt$p.value,
               note = paste(
                   "secondary relative contrast only; quartile boundaries are",
                   "not biological or absolute-PVE cutoffs"))
} else {
    data.table(model = "top_bottom_quartiles_SECONDARY", note = "too few loci")
}
write_atomic(secondary, file.path(res_dir, "burden-extreme-groups.tsv"))

## Emit the burden table under the names `config/meqtl_parameters.yml:
## vmr_aggregation_metrics` declares, which is what 04_check_burden.R reads.
## The internal short names are kept for the model formulae above and dropped
## here, so the written table has exactly one name per quantity.
data.table::setnames(burden,
                     c("n_sig", "prop_sig"),
                     c("n_cpgs_with_sig_meqtl", "proportion_cpgs_with_sig_meqtl"))
write_atomic(burden, file.path(res_dir, "vmr-meqtl-burden.tsv"))

writeLines(c(
    "Interpretation constraints carried by this table:",
    "  - meQTL mapping used the SAME donors as the source estimates. This is",
    "    convergent evidence, NOT independent replication.",
    "  - Burden denominators are TESTED CpGs; see excluded-cpgs.tsv and",
    "    vmrs-without-tested-cpgs.tsv for what is not in them.",
    "  - Public positive-only meQTL resources cannot provide an external",
    "    gradient: absence from a positive list is not a tested negative.",
    sprintf("  - Genomic inflation lambda = %.3f.", lambda)
), file.path(res_dir, "interpretation-constraints.txt"))

print(primary[term == "local_snp_contribution_score_z"])
