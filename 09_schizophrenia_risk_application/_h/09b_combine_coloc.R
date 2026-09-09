#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- pool the per-chromosome coloc runs ####
##
## Usage:
##   Rscript _h/09b_combine_coloc.R --run-id scz-AA-caudate-YYYYMMDD
##
## Reconciles the coloc fan-out and rolls the region-level posteriors up to the
## locus level. A chromosome that produced no evaluable region is a recordable
## outcome; a chromosome that produced no FILE is an unexplained computational
## failure and stops the run (AGENTS.md 9).
##
## Colocalization posteriors are NOT multiplicity-corrected. PP4 is already a
## posterior probability under stated priors, not a p-value, and applying FDR to
## it would be a category error. The prespecified decision rule is the PP4 and
## PP4/PP3 thresholds in config, applied per region.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

MODULE <- "09_schizophrenia_risk_application"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mf <- function(field) {
    v <- manifest[["value"]][manifest[["field"]] == field]
    if (length(v) == 0) NA_character_ else v[1]
}
smoke <- identical(mf("smoke_run"), "TRUE")

in_dir <- file.path(run_dir, "coloc", "abf")
expected <- paste0("chr", 1:22)
present <- sub("\\.tsv\\.gz$", "", basename(
    list.files(in_dir, pattern = "^chr[0-9]+\\.tsv\\.gz$")))
reconcile(expected = expected, completed = present,
          failed = setdiff(expected, present),
          run = list(dir = run_dir), allow_failures = smoke)

parts <- lapply(file.path(in_dir, paste0(present, ".tsv.gz")), function(f) {
    dt <- fread(f)
    if (nrow(dt) == 0) NULL else dt
})
res <- rbindlist(Filter(Negate(is.null), parts), use.names = TRUE, fill = TRUE)

if (nrow(res) == 0) {
    write_atomic(data.table(), file.path(run_dir, "results", "coloc-abf.tsv"))
    write_atomic(data.table(), file.path(run_dir, "results", "coloc-locus-summary.tsv"))
    message("[09] no coloc region was prepared on any autosome")
    quit(status = 0)
}

fwrite(res, file.path(run_dir, "results", "coloc-abf.tsv.gz"),
       sep = "\t", compress = "gzip")

## Locus level, split by arm: the gate and the prioritization rule both read
## the gate-eligible arms only, and keeping the cross-ancestry arm in its own
## columns is what stops it leaking into a claim.
locus <- res[, .(
    n_regions = .N,
    n_regions_evaluated = sum(status == "OK"),
    n_regions_unevaluable = sum(status != "OK"),
    max_pp4 = if (any(status == "OK")) max(PP4, na.rm = TRUE) else NA_real_,
    n_colocalized = sum(colocalized, na.rm = TRUE),
    n_claimable = sum(claimable_colocalization, na.rm = TRUE),
    top_phenotype = if (any(status == "OK"))
        phenotype_id[which.max(fifelse(status == "OK", PP4, -Inf))] else NA_character_,
    top_context = if (any(status == "OK"))
        qtl_context[which.max(fifelse(status == "OK", PP4, -Inf))] else NA_character_
), by = .(locus_id, arm, gate_eligible, ld_ancestry_matched, arm_status)]
setorder(locus, -max_pp4, locus_id, arm)
write_atomic(locus, file.path(run_dir, "results", "coloc-locus-summary.tsv"))

## One row per locus for the prioritization join and the gate: gate-eligible
## arms only.
gated <- res[gate_eligible == TRUE & ld_ancestry_matched == TRUE]
per_locus <- if (nrow(gated) == 0) {
    data.table(locus_id = unique(res$locus_id),
               max_pp4_gate_eligible = NA_real_,
               has_claimable_colocalization = FALSE)
} else {
    gated[, .(max_pp4_gate_eligible =
                  if (any(status == "OK")) max(PP4, na.rm = TRUE) else NA_real_,
              has_claimable_colocalization = any(claimable_colocalization, na.rm = TRUE)),
          by = locus_id]
}
write_atomic(per_locus, file.path(run_dir, "results", "coloc-locus-gate-eligible.tsv"))

message("[09] coloc: ", nrow(res), " regions over ", uniqueN(res$locus_id),
        " loci; ", sum(res$claimable_colocalization, na.rm = TRUE),
        " claimable colocalizations")
