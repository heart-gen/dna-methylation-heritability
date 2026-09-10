#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- link PGC3 loci to corrected VMRs ####
##
## Usage:
##   Rscript _h/02_link_loci_to_vmrs.R --run-id scz-AA-caudate-YYYYMMDD
##
## The linkage is positional and prespecified: a VMR is reachable from a locus
## if it lies inside the locus interval widened by locus_flank_bp. The flank is
## config/meqtl_parameters.yml:cis_window_bp, so a risk variant can only reach a
## CpG that Module 05 actually tested it against -- widening it here would
## nominate pairs the meQTL scan never evaluated.
##
## This stage also freezes the tested universe. AGENTS.md 7.8 requires the
## risk-variant x CpG tests to form their own FDR family, and a family whose
## membership is decided after the p-values are seen is not a family.

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

loci <- fread(file.path(run_dir, "results", "scz-loci.tsv"))
if (nrow(loci) == 0) stop("No loci; run 01_define_scz_loci.R first")

## ------------------------------------------------- the corrected VMR universe
## Module 02's table is the VMR universe for every downstream biological
## module: it carries the boundaries, the relative score, and the eligibility
## flag in one place, and its loader refuses superseded estimator columns.
lgc <- load_local_genetic_control(mf("upstream_local_genetic_variance_run_id"),
                                  region = region, cohort = cohort,
                                  eligible_only = TRUE)
vmrs <- as.data.table(lgc)[, .(vmr_id, chrom, start, end, n_cpgs, n_variants,
                               vmr_set_id,
                               local_snp_contribution_score,
                               local_snp_contribution_score_z,
                               local_snp_contribution_quartile)]
vmrs[, chrom := paste0("chr", sub("^chr", "", chrom))]

win <- loci[, .(chrom, start = window_start, end = window_end, locus_id,
                locus_start = start, locus_end = end,
                pgc3_index_pvalue, index_snp)]
setkey(win, chrom, start, end)

links <- foverlaps(vmrs[, .(vmr_id, chrom, start, end, n_cpgs, n_variants,
                            local_snp_contribution_score,
                            local_snp_contribution_score_z,
                            local_snp_contribution_quartile)],
                   win, by.x = c("chrom", "start", "end"),
                   type = "any", nomatch = NULL)
setnames(links, c("i.start", "i.end"), c("vmr_start", "vmr_end"))
links[, `:=`(window_start = start, window_end = end, start = NULL, end = NULL)]
## Distance from the VMR to the published interval, 0 when it lies inside. Kept
## so a reader can see which links are inside the locus and which are in the
## cis flank.
links[, distance_to_locus := pmax(0L,
    pmax(locus_start - vmr_end, vmr_start - locus_end))]
setorder(links, locus_id, distance_to_locus, vmr_id)
write_atomic(links, file.path(run_dir, "results", "locus-vmr-links.tsv"))

## ------------------------------------------------------------ tested universe
universe <- data.table(
    run_id = opts$run_id, cohort = cohort, region = region,
    vmr_set_id = if (nrow(vmrs)) vmrs$vmr_set_id[1] else NA_character_,
    n_loci_published = nrow(loci),
    n_loci_with_linked_vmr = uniqueN(links$locus_id),
    n_vmrs_total = nrow(vmrs),
    n_vmrs_linked = uniqueN(links$vmr_id),
    n_locus_vmr_links = nrow(links),
    locus_flank_bp = as.integer(scz$locus_definition$locus_flank_bp),
    smoke_run = smoke
)
write_atomic(universe, file.path(run_dir, "results", "tested-universe.tsv"))

## The background for every VMR-level enrichment: eligible VMRs NOT reachable
## from any schizophrenia locus. Written out so the comparison set is a
## recorded object rather than a filter re-derived in three later scripts.
background <- vmrs[!vmr_id %in% links$vmr_id]
background[, scz_linked := FALSE]
linked <- vmrs[vmr_id %in% links$vmr_id]
linked[, scz_linked := TRUE]
write_atomic(rbind(linked, background),
             file.path(run_dir, "results", "vmr-scz-linkage.tsv"))

print(universe[, .(n_loci_with_linked_vmr, n_vmrs_linked, n_locus_vmr_links)])
message("[09] linked ", uniqueN(links$vmr_id), " VMRs to ",
        uniqueN(links$locus_id), " loci")
