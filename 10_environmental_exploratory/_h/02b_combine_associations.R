#!/usr/bin/env Rscript
#### 10_environmental_exploratory -- pool chromosomes and correct once ####
##
## Usage:
##   Rscript _h/02b_combine_associations.R --run-id env-AA-caudate-YYYYMMDD
##
## AGENTS.md 10.3: "FDR families are not combined after inspection." One BH
## family per exposure, over the VMRs tested for that exposure, applied here and
## nowhere else. Correcting inside the per-chromosome tasks would give 22
## families per exposure; correcting across exposures afterwards would merge
## families that were declared separate.
##
## AGENTS.md 9: a SLURM array that exits 0 on every task is not evidence that
## every task produced output, so the 22 chromosomes are reconciled explicitly.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv(
    "V2_RUN_CODE",
    file.path(Sys.getenv("V2_REPO_ROOT", "."), "10_environmental_exploratory", "_h")),
    "run_config.R"))

suppressPackageStartupMessages({
    library(data.table)
})

MODULE <- "10_environmental_exploratory"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
env <- load_run_config("environmental", run_dir)
alpha <- as.numeric(env$testing$fdr_alpha)

chrom_dir <- file.path(run_dir, "results", "per-chrom")
expected <- as.character(1:22)
present <- sub("^assoc-chr", "",
               sub("\\.tsv$", "",
                   list.files(chrom_dir, pattern = "^assoc-chr.*\\.tsv$")))
parts <- lapply(present, function(cc) {
    f <- file.path(chrom_dir, paste0("assoc-chr", cc, ".tsv"))
    if (file.size(f) <= 1) return(NULL)
    fread(f)
})
names(parts) <- present
empty <- present[vapply(parts, function(x) is.null(x) || !nrow(x), logical(1))]

reconcile(expected = expected,
          completed = setdiff(present, empty),
          excluded = empty,
          run = list(dir = run_dir))

res <- rbindlist(Filter(Negate(is.null), parts), fill = TRUE)
if (!nrow(res)) stop("No association rows in any chromosome")

ok <- res[status == "ok"]
if (!nrow(ok)) stop("No VMR x exposure model fitted successfully")

## One row per VMR x exposure carries the joint test; the per-level coefficient
## rows are retained beside it but are not a second family.
per_vmr <- unique(ok[, .(cohort, region, run_id, vmr_id, chrom, start, end,
                         n_cpgs, exposure, stratum, p_joint, df_num, n_used)])
if (anyDuplicated(per_vmr[, .(vmr_id, exposure, stratum)])) {
    stop("A VMR x exposure x stratum triple produced more than one joint p-value")
}

## BH within exposure AND stratum. The same exposure tested pooled and within
## cases is two questions on two donor sets, not one family split in half, so
## `by = .(exposure, stratum)` is the whole point of this line.
per_vmr[, fdr := stats::p.adjust(p_joint, method = "BH"),
        by = .(exposure, stratum)]
per_vmr[, significant := fdr <= alpha]

fam <- per_vmr[, .(n_tested_vmrs = .N,
                   n_significant = sum(significant, na.rm = TRUE),
                   min_p = min(p_joint, na.rm = TRUE),
                   median_n_donors = as.integer(stats::median(n_used))),
               by = .(exposure, stratum)]
fam[, `:=`(fdr_alpha = alpha, fdr_method = "BH",
           cohort = mf("cohort"), region = mf("region"), run_id = opts$run_id)]

## Every emitted row states what it may be used for, so a downstream consumer
## inherits the constraint from the data rather than from prose (AGENTS.md 2.3).
per_vmr[, `:=`(exploratory_supplement_only = TRUE,
               causal_interpretation_allowed = FALSE)]
fam[, `:=`(exploratory_supplement_only = TRUE,
           causal_interpretation_allowed = FALSE)]

write_atomic(per_vmr, file.path(run_dir, "results",
                                "vmr-exposure-association.tsv"))
write_atomic(ok, file.path(run_dir, "results",
                           "vmr-exposure-association-terms.tsv"))
write_atomic(fam, file.path(run_dir, "results", "fdr-families.tsv"))

n_families <- nrow(fam)
declared <- strsplit(mf("eligible_exposures"), ",", fixed = TRUE)[[1]]
declared <- declared[nzchar(declared)]
if (n_families != length(declared)) {
    stop("Built ", n_families, " BH families for ", length(declared),
         " declared eligible exposure-stratum pairs. The family structure must ",
         "match what stage 01 prespecified (AGENTS.md 10.3).")
}

append_manifest(list(dir = run_dir), list(
    n_fdr_families = as.character(n_families),
    n_tested_vmr_exposure_pairs = as.character(nrow(per_vmr)),
    n_significant_pairs = as.character(sum(per_vmr$significant, na.rm = TRUE)),
    n_failed_models = as.character(nrow(res[status %in% c("model_failed",
                                                          "phenotype_unavailable")]))
))
print(fam[, .(exposure, stratum, n_tested_vmrs, n_significant, min_p)])
