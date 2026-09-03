#!/usr/bin/env Rscript
##
## Post-hoc QC: where does the v1 enrichment signal go under the v2 adjustment?
##
##   Rscript _h/06_qc_adjustment_ladder.R --run-id rra-AA-caudate-20260902
##   (no --run-id: all three regions of the current production date)
##
## v1 (meqtl-validation/02_vmr_meqtl_burden) reported H3K9me3 enrichment at
## OR 1.36 in caudate and hippocampus and quiescent-chromatin enrichment in all
## three regions at FDR <= 2e-10, adjusting for VMR length and local SNP count
## only. v2 returns nulls for both. The two analyses differ in several ways, but
## the adjustment set is the one that can be isolated on the sealed feature
## table, so this script does that: refit the primary model along a nested
## ladder from the v1-like adjustment up to the full v2 set, holding the
## predictor, outcome scale, family and robust SEs fixed.
##
## This is post-hoc analysis OF a sealed run, not part of it. It writes to
## _m/qc/{run_id}/, which is gitignored and regenerable.
##
## Read the output as a decomposition, NOT as a menu to pick from. A smaller
## adjustment set is not automatically the better estimate; the question the
## ladder answers is WHICH covariates carry the attenuation, so the confounder
## vs mediator/collider judgement can be made explicitly per covariate.

suppressPackageStartupMessages({
    library(data.table); library(sandwich); library(lmtest)
})
source(file.path(Sys.getenv("V2_REPO", "."), "00_shared", "load.R"))

MODULE <- "04_repeat_repressive_architecture"
PRED   <- "local_snp_contribution_score_z"

## The ladder. Each rung ADDS to the one above it, so the change in the estimate
## between adjacent rungs is attributable to the covariates newly introduced.
## Terms match _h/02_test_association.R exactly, including log(vmr_length) --
## fitting raw length here would not reproduce the sealed run.
LADDER <- list(
    v1_like   = c("log(vmr_length)", "tested_snp_count"),
    plus_cpg  = c("cpg_count", "cpg_density"),
    plus_seq  = c("gc_content"),
    plus_meth = c("mean_methylation", "methylation_variance", "wgbs_coverage"),
    plus_tech = c("snp_proximal_frac", "mappability", "segdup_frac",
                  "problematic_frac"),
    full_v2   = c("factor(broad_genomic_annotation)", "cell_composition_r2")
)
OUTCOMES <- c("h3k9me3_frac", "quiescent_frac", "line_l1_frac")

fit_rung <- function(d, outcome, terms) {
    form <- as.formula(paste(outcome, "~", PRED, "+",
                             paste(terms, collapse = " + ")))
    fit <- tryCatch(glm(form, data = d, family = quasibinomial(link = "logit")),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    ct <- tryCatch(coeftest(fit, vcov. = vcovHC(fit, type = "HC3")),
                   error = function(e) NULL)
    if (is.null(ct) || !PRED %in% rownames(ct)) return(NULL)
    data.table(n_fitted = nobs(fit), estimate = ct[PRED, 1],
               se = ct[PRED, 2], p = ct[PRED, 4])
}

run_one <- function(run_id) {
    d <- fread(file.path(repo_root(), MODULE, "_m", "runs", run_id,
                         "results", "vmr-features.tsv"))
    region <- sub("^rra-[A-Za-z]+-([a-z]+)-.*$", "\\1", run_id)
    rbindlist(lapply(OUTCOMES, function(o) {
        cum <- character(0)
        rbindlist(lapply(names(LADDER), function(nm) {
            cum <<- c(cum, LADDER[[nm]])
            r <- fit_rung(d, o, cum)
            if (is.null(r)) return(NULL)
            cbind(data.table(run_id = run_id, region = region, outcome = o,
                             adjustment = nm,
                             added = paste(LADDER[[nm]], collapse = ", ")), r)
        }))
    }))
}

opts <- commandArgs(trailingOnly = TRUE)
run_ids <- if (length(opts) >= 2 && opts[1] == "--run-id") opts[2] else {
    dirs <- list.dirs(file.path(repo_root(), MODULE, "_m", "runs"),
                      recursive = FALSE, full.names = FALSE)
    grep("^rra-AA-[a-z]+-[0-9]{8}$", dirs, value = TRUE)
}
if (length(run_ids) == 0) stop("no run ids resolved")

res <- rbindlist(lapply(run_ids, run_one))

## Order matters here: the rungs are cumulative, so a table sorted
## alphabetically reads as noise rather than as an attenuation sequence.
res[, adjustment := factor(adjustment, levels = names(LADDER))]
setorder(res, run_id, outcome, adjustment)

## Attenuation relative to the v1-like rung, which is the comparison the
## question is about. Reported as a retained fraction of the v1-like estimate.
res[, retained_vs_v1 := estimate / estimate[adjustment == "v1_like"],
    by = .(run_id, outcome)]

qc_dir <- file.path(repo_root(), MODULE, "_m", "qc", run_ids[1])
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(res, file.path(qc_dir, "adjustment-ladder.tsv"), sep = "\t")

for (o in OUTCOMES) {
    cat("\n===", o, "\n")
    print(dcast(res[outcome == o], adjustment ~ region,
                value.var = c("estimate", "p")), digits = 3)
}
cat("\nwrote", file.path(qc_dir, "adjustment-ladder.tsv"), "\n")
