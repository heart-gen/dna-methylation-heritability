#!/usr/bin/env Rscript
#### 05_cpg_meqtl_burden -- pool chromosomes and apply the region FDR ####
##
## Usage:
##   Rscript _h/02b_combine_meqtl.R --run-id cmb-AA-caudate-20260823
##
## config/meqtl_parameters.yml sets `fdr_family: per_brain_region`. The FDR
## family therefore spans every autosome of one region, and the correction can
## only be applied AFTER the per-chromosome array has finished. Applying Storey
## inside each array task would create 22 unrelated families and make the burden
## fraction depend on how the work happened to be split across jobs.
##
## This stage also closes the tested-CpG accounting. AGENTS.md 7.5 requires
## tested CpGs to be reported separately from prepared-but-untested ones, and
## every concordance denominator to be audited. A CpG that was prepared but had
## no variant in its cis window is untestable, not negative; counting it as a
## CpG without meQTL support would understate the burden of exactly those VMRs
## that sit in variant-poor regions -- which are not a random subset.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

MODULE <- "05_cpg_meqtl_burden"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]; if (length(v) == 0) NA_character_ else v[1]
}
smoke <- identical(mval("smoke_run"), "TRUE")
meqtl <- load_config("meqtl_parameters")

res_dir <- file.path(run_dir, "results")
map_dir <- file.path(res_dir, "meqtl")
if (!dir.exists(map_dir)) stop("No mapping output under ", map_dir)

## --------------------------------------------------------- reconcile chroms
## The array covers the autosomes the config declares, and each one must have
## either results or an explicit skip marker.
expected_chroms <- paste0("chr", config_get(load_config("thresholds"),
                                            "chromosomes.primary"))
## A smoke run that stage 01 restricted to a subset of autosomes reconciles
## against that recorded subset, never against the full policy list. Honoured
## only when the manifest also says `smoke_run = TRUE`.
if (smoke && !is.na(mval("smoke_chroms")) && nzchar(mval("smoke_chroms"))) {
    expected_chroms <- trimws(strsplit(mval("smoke_chroms"), ",")[[1]])
    message("[05] SMOKE: reconciling against ",
            paste(expected_chroms, collapse = ","), " only")
}
res_files <- list.files(map_dir, pattern = "\\.cis_qtl\\.tsv\\.gz$", full.names = TRUE)
have <- sub("\\.cis_qtl\\.tsv\\.gz$", "", basename(res_files))
skipped <- sub("\\.skipped$", "", list.files(map_dir, pattern = "\\.skipped$"))
reconcile(expected = expected_chroms, completed = have, qc_failed = skipped,
          run = list(dir = run_dir), allow_failures = smoke)

if (length(res_files) == 0) stop("Every autosome was skipped; nothing to combine")

cpg_res <- rbindlist(lapply(res_files, function(f) {
    d <- fread(f)
    d[, chrom := sub("\\.cis_qtl\\.tsv\\.gz$", "", basename(f))]
    d
}), fill = TRUE)

## tensorqtl names the beta-approximated permutation p-value `pval_beta` and
## falls back to `pval_perm` when the beta fit did not converge. Prefer the
## beta p-value where it exists, and RECORD which one each CpG used rather than
## silently mixing two scales in one column.
cpg_res[, pval_source := fifelse(!is.na(pval_beta), "pval_beta", "pval_perm")]
cpg_res[, pvalue := fifelse(!is.na(pval_beta), pval_beta, pval_perm)]
cpg_res <- cpg_res[is.finite(pvalue)]
if (nrow(cpg_res) == 0) stop("No CpG produced a finite permutation p-value")

## ------------------------------------------------------------- region FDR
## Storey q-value where qvalue is available; Benjamini-Hochberg otherwise, with
## the method actually used recorded in the output. BH is conservative relative
## to Storey, so a fallback can only understate the burden, never inflate it.
##
## Three sources of Storey q-values are tried in order:
##   1. the Bioconductor `qvalue` package, if this R environment has it;
##   2. `_h/storey_qvalue.py`, which uses py_qvalue under the `genomics`
##      environment -- `epigenomics` (this stage) has neither package, and the
##      module already crosses into `genomics` for tensorqtl, so shelling out
##      is cheaper than adding a package to a shared environment;
##   3. BH.
##
## The distinction matters here rather than being a formality: on the AA caudate
## chr22 smoke, pi0 = 0.217 and Storey called 2,098 CpGs at FDR 0.05 against
## BH's 1,636 -- a 28% difference in the module's denominator-normalised
## endpoint.
storey_via_python <- function(p) {
    py <- "/projects/p32505/opt/envs/genomics/bin/python"
    script <- file.path(repo_root(), "05_cpg_meqtl_burden", "_h",
                        "storey_qvalue.py")
    if (!file.exists(py) || !file.exists(script)) return(NULL)
    out <- suppressWarnings(tryCatch(
        system2(py, shQuote(script), stdout = TRUE, stderr = FALSE,
                input = format(p, scientific = TRUE, digits = 17)),
        error = function(e) NULL))
    if (is.null(out) || !is.null(attr(out, "status")) ||
        length(out) != length(p) + 1L || !startsWith(out[1], "pi0\t")) {
        return(NULL)
    }
    q <- as.numeric(out[-1])
    if (anyNA(q)) return(NULL)
    list(q = q, pi0 = as.numeric(sub("^pi0\t", "", out[1])))
}

fdr_method <- meqtl$mapping$fdr_method
used <- NULL
if (identical(fdr_method, "storey_qvalue")) {
    if (requireNamespace("qvalue", quietly = TRUE)) {
        cpg_res[, qvalue := qvalue::qvalue(pvalue)$qvalues]
        used <- "storey_qvalue_bioconductor"
    } else {
        st <- storey_via_python(cpg_res$pvalue)
        if (!is.null(st)) {
            cpg_res[, qvalue := st$q]
            message("[05] Storey q-values via py_qvalue; pi0 = ",
                    signif(st$pi0, 4))
            used <- "storey_qvalue_py_qvalue"
        }
    }
}
if (is.null(used)) {
    cpg_res[, qvalue := stats::p.adjust(pvalue, method = "BH")]
    used <- "benjamini_hochberg"
    if (identical(fdr_method, "storey_qvalue")) {
        warning("config requests storey_qvalue but neither the qvalue package ",
                "nor py_qvalue could be used; fell back to BH, which is ",
                "conservative. Recorded in the output as fdr_method_used.",
                call. = FALSE)
    }
}
cpg_res[, `:=`(fdr_family = "per_brain_region", fdr_method_used = used,
               region = mval("region"), population = mval("cohort"))]
write_atomic(cpg_res, file.path(res_dir, "cpg-meqtl-results.tsv"))

## ---------------------------------------------------- denominator accounting
untested <- rbindlist(lapply(
    list.files(map_dir, pattern = "\\.untested-cpgs\\.tsv$", full.names = TRUE),
    fread), fill = TRUE)
if (nrow(untested)) write_atomic(untested, file.path(res_dir, "untested-cpgs.tsv"))

member <- fread(file.path(res_dir, "tested-cpg-membership.tsv"))
prepared <- nrow(member)
tested <- nrow(cpg_res)
audit <- data.table(
    region = mval("region"), population = mval("cohort"),
    n_prepared_cpgs = prepared,
    n_tested_cpgs = tested,
    n_prepared_but_untested = nrow(untested),
    n_unaccounted = prepared - tested - nrow(untested),
    n_significant_cpgs = sum(cpg_res$qvalue <= meqtl$mapping$fdr_threshold),
    fdr_threshold = meqtl$mapping$fdr_threshold,
    fdr_method_used = used,
    chromosomes_mapped = length(res_files),
    chromosomes_skipped = length(skipped)
)
write_atomic(audit, file.path(res_dir, "meqtl-denominator-audit.tsv"))
print(audit)

## Prepared must equal tested plus untested. A nonzero remainder means CpGs
## vanished between preparation and mapping, which invalidates every burden
## fraction downstream, so it stops the run rather than being noted.
if (audit$n_unaccounted != 0L && !smoke) {
    stop("CpG accounting does not balance: ", prepared, " prepared, ", tested,
         " tested, ", nrow(untested), " explicitly untested, ",
         audit$n_unaccounted, " unaccounted.")
}
message("[05] combined ", length(res_files), " chromosome(s); ",
        audit$n_significant_cpgs, " CpGs with meQTL support at FDR ",
        meqtl$mapping$fdr_threshold)
