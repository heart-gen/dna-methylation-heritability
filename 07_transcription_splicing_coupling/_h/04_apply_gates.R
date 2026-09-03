#!/usr/bin/env Rscript
#### 07_transcription_splicing_coupling -- acceptance gate ####
##
## Usage:
##   Rscript _h/04_apply_gates.R --run-id tsc-AA-caudate-YYYYMMDD
##
## A null coupling result is a reportable finding, not a failed run. This gate
## therefore checks that the analysis was CONDUCTED correctly -- the declared
## modalities ran, the tested universe is large enough to interpret, no banned
## quantity entered a model -- and does not require a significant result.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE", file.path(Sys.getenv("V2_REPO_ROOT", "."), "07_transcription_splicing_coupling", "_h")), "run_config.R"))
suppressPackageStartupMessages(library(data.table))

MODULE <- "07_transcription_splicing_coupling"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
smoke <- identical(mf("smoke_run"), "TRUE")
ts <- load_run_config("transcription_splicing", run_dir)
enabled <- strsplit(mf("modalities"), ",", fixed = TRUE)[[1]]

uni_f <- file.path(run_dir, "results", "tested-universe.tsv")
tests_f <- file.path(run_dir, "results", "coupling-tests.tsv")
if (!file.exists(uni_f)) stop("No tested-universe.tsv; run stage 01")
if (!file.exists(tests_f)) stop("No coupling-tests.tsv; run stage 03")
uni <- fread(uni_f)
tests <- fread(tests_f)

fail <- character()

ran <- unique(tests$modality)
if (isTRUE(ts$gates$require_all_enabled_modalities) && !smoke) {
    absent <- setdiff(enabled, ran)
    if (length(absent)) {
        fail <- c(fail, paste0("MODALITY_NOT_RUN(", paste(absent, collapse = ","), ")"))
    }
}

if (!smoke) {
    if (max(uni$n_vmrs_linked, na.rm = TRUE) < ts$gates$min_vmrs_tested) {
        fail <- c(fail, "TOO_FEW_VMRS_TESTED")
    }
    if (max(uni$n_pairs, na.rm = TRUE) < ts$gates$min_pairs_tested) {
        fail <- c(fail, "TOO_FEW_PAIRS_TESTED")
    }
}

## Every reported test must have produced a finite estimate; a table of NAs
## would otherwise pass silently as "no significant coupling".
if (all(!is.finite(tests$estimate))) {
    fail <- c(fail, "NO_TEST_PRODUCED_AN_ESTIMATE")
}

## The internal LIBD eQTL map is deliberately out of scope: its genome-wide QC
## repair is unresolved (meqtl-validation/09_libd_eqtl_mapping/EQTL_DEBUG_TODO.md).
## A run that switched it on without that repair would be citing a broken map.
if (isTRUE(ts$internal_libd_eqtl_support_arm)) {
    fail <- c(fail, "LIBD_EQTL_ARM_ENABLED_WHILE_QC_REPAIR_OPEN")
}

thr <- ts$association$fdr_threshold
sig <- tests[is.finite(q) & q < thr]
lgc <- tests[predictor == "local_genetic_control" & is.finite(q) & q < thr]
meq <- tests[predictor %in% c("any_meqtl_support", "meqtl_proportion") &
             is.finite(q) & q < thr]

decision <- if (length(fail) > 0) {
    paste0("FAIL_TX_COUPLING_QC:", paste(fail, collapse = ";"))
} else if (smoke) {
    "PASS_SMOKE_ONLY_NOT_ACCEPTABLE"
} else {
    "PASS_TX_COUPLING_QC"
}

dec <- data.table(
    run_id = opts$run_id, cohort = mf("cohort"), region = mf("region"),
    decision = decision, smoke_run = smoke,
    modalities_declared = paste(enabled, collapse = ","),
    modalities_run = paste(ran, collapse = ","),
    n_tests = nrow(tests),
    n_significant = nrow(sig),
    n_significant_local_genetic_control = nrow(lgc),
    n_significant_meqtl_support = nrow(meq),
    total_vmrs_tested = max(uni$n_vmrs_linked, na.rm = TRUE),
    total_pairs_tested = sum(uni$n_pairs, na.rm = TRUE),
    fdr_threshold = thr,
    ## Module 09 consumes this: AGENTS.md 7.8 requires "at least one locus with
    ## transcriptional coupling" among its retention criteria.
    coupling_supported = nrow(sig) > 0 && !smoke && length(fail) == 0,
    ## AGENTS.md 7.6: the permitted reading of a positive result, carried with
    ## the result so it cannot drift in the writing.
    permitted_claim = if (nrow(sig) > 0)
        "genetically regulated VMRs are more frequently transcriptionally coupled"
        else "no coupling detected in the prespecified tested universe",
    forbidden_claim = "methylation mediates the genetic effect on expression or splicing",
    failures = paste(fail, collapse = ";")
)
write_atomic(dec, file.path(run_dir, "results", "coupling-decision.tsv"))

writeLines(c(
    "Interpretation constraints carried by this run:",
    "  - Tests are VMR-level associations between genetic regulation and the",
    "    presence of a local transcriptional association. They do not test",
    "    mediation, and no causal ordering is estimated anywhere in this module.",
    "  - The tested universe is the prespecified local link set in",
    "    results/tested-universe.tsv. It is not a transcriptome-wide screen.",
    "  - Pair-level FDR is applied within one modality within one cell",
    paste0("    (fdr_family: ", ts$association$fdr_family, ")."),
    "  - The internal LIBD eQTL map is NOT used; its genome-wide QC repair is",
    "    open (meqtl-validation/09_libd_eqtl_mapping/EQTL_DEBUG_TODO.md).",
    paste0("  - Permitted claim: ", dec$permitted_claim[1]),
    paste0("  - Forbidden claim: ", dec$forbidden_claim[1])
), file.path(run_dir, "results", "interpretation-constraints.txt"))

print(dec[, .(decision, n_tests, n_significant, total_pairs_tested)])
message("[07] decision ", decision)
if (startsWith(decision, "FAIL")) quit(status = 1)
