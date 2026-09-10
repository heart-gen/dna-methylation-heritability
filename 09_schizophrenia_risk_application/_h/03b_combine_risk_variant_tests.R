#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- pool chromosomes, correct once ####
##
## Usage:
##   Rscript _h/03b_combine_risk_variant_tests.R --run-id scz-AA-caudate-YYYYMMDD
##
## AGENTS.md 7.8 requires the risk-variant x CpG tests to keep their own FDR
## family, and config/schizophrenia.yml names it
## `scz_risk_variant_cpg_pairs_per_region`. That family spans all autosomes of
## one region, so the correction has to happen after the per-chromosome tasks
## are pooled -- exactly the structure 05_cpg_meqtl_burden uses, and for the
## same reason: an FDR family whose membership depends on how the work was
## parallelised is not a family.
##
## This is also where the run reconciles its fan-out. A chromosome that produced
## no pairs is a recordable outcome; a chromosome that produced no FILE is an
## unexplained computational failure and stops the run (AGENTS.md 9).

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "09_schizophrenia_risk_application", "_h")),
    "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "09_schizophrenia_risk_application"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mf("cohort"); region <- mf("region")
smoke <- identical(mf("smoke_run"), "TRUE")

scz <- load_run_config("schizophrenia", run_dir)
alpha <- as.numeric(scz$testing$fdr_alpha)
method <- scz$testing$fdr_method

in_dir <- file.path(run_dir, "results", "risk_variant_tests")
expected <- paste0("chr", 1:22)
present <- sub("\\.tsv\\.gz$", "", basename(
    list.files(in_dir, pattern = "^chr[0-9]+\\.tsv\\.gz$")))
reconcile(expected = expected, completed = present,
          failed = setdiff(expected, present),
          run = list(dir = run_dir), allow_failures = smoke)

parts <- lapply(file.path(in_dir, paste0(present, ".tsv.gz")), function(f) {
    dt <- fread(f)
    if (nrow(dt) == 0) return(NULL)
    dt
})
pairs <- rbindlist(Filter(Negate(is.null), parts), use.names = TRUE, fill = TRUE)

if (nrow(pairs) == 0) {
    ## No pair anywhere is a reportable null, not a crash. Write the empty
    ## tables the later stages expect so they fail on their own gates rather
    ## than on a missing file.
    empty <- data.table()
    write_atomic(empty, file.path(run_dir, "results", "risk-variant-cpg-tests.tsv"))
    write_atomic(empty, file.path(run_dir, "results", "locus-meqtl-support.tsv"))
    stop("No risk-variant x CpG pair was tested on any autosome. Check the ",
         "Module 05 nominal outputs and the locus-VMR links before rerunning.")
}

## One correction, over the whole family, once.
pairs[, qvalue := p.adjust(pval_nominal, method = method)]
pairs[, significant := is.finite(qvalue) & qvalue < alpha]
pairs[, fdr_family := scz$testing$fdr_family]
pairs[, fdr_method := method]
setorder(pairs, qvalue, pval_nominal)
fwrite(pairs, file.path(run_dir, "results", "risk-variant-cpg-tests.tsv.gz"),
       sep = "\t", compress = "gzip")

## Locus-level support: the unit the retention criterion and the prioritization
## rule are both stated in.
locus <- pairs[, .(
    n_risk_variants_tested = uniqueN(risk_variant_id),
    n_cpgs_tested = uniqueN(cpg_id),
    n_vmrs_tested = uniqueN(vmr_id),
    n_pairs_tested = .N,
    n_pairs_significant = sum(significant),
    n_vmrs_significant = uniqueN(vmr_id[significant]),
    min_pvalue = min(pval_nominal, na.rm = TRUE),
    min_qvalue = min(qvalue, na.rm = TRUE),
    max_abs_slope = max(abs(slope), na.rm = TRUE),
    has_significant_risk_variant_cpg_meqtl = any(significant)
), by = .(locus_id, index_snp, chrom)]
setorder(locus, min_qvalue, locus_id)
write_atomic(locus, file.path(run_dir, "results", "locus-meqtl-support.tsv"))

summary_dt <- data.table(
    run_id = opts$run_id, cohort = cohort, region = region,
    fdr_family = scz$testing$fdr_family, fdr_method = method, fdr_alpha = alpha,
    n_chromosomes_with_pairs = uniqueN(pairs$chrom),
    n_loci_tested = nrow(locus),
    n_loci_with_significant_meqtl = sum(locus$has_significant_risk_variant_cpg_meqtl),
    n_pairs_tested = nrow(pairs),
    n_pairs_significant = sum(pairs$significant),
    n_vmrs_with_significant_meqtl = uniqueN(pairs$vmr_id[pairs$significant]),
    smoke_run = smoke
)
write_atomic(summary_dt, file.path(run_dir, "results", "risk-variant-test-summary.tsv"))
print(summary_dt[, .(n_loci_tested, n_loci_with_significant_meqtl,
                     n_pairs_tested, n_pairs_significant)])
message("[09] ", sum(locus$has_significant_risk_variant_cpg_meqtl), "/",
        nrow(locus), " loci carry CpG meQTL support at FDR < ", alpha)
