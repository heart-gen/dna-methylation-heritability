#!/usr/bin/env Rscript
#### 09_schizophrenia_gwas_loci -- open the targeted EA check ####
##
## Usage:
##   Rscript _h/19_ea_targeted_new_run.R [--allow-unlocked] [--run-id ID]
##
## Opens `scz-all_individuals.EA-crossregion-{date}` and writes the targets it
## will test, before any EA genotype or methylation is read:
##
##   results/targets.tsv   one row per (region, locus, risk variant, CpG): every
##                         pair significant in the AA primary at that region's
##                         prespecified illustrative loci (prioritized == TRUE),
##                         with the AA slope, SE and p it is checked against
##   results/loci.tsv      the loci, their lead pair (minimum AA p)
##
## Everything is read from the accepted runs named in config/scz_ea_targeted.yml,
## and each is checked against the module README's accepted-runs table, so a
## superseded run cannot be the source of the targets.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

MODULE <- "09_schizophrenia_gwas_loci"
COHORT <- "all_individuals.EA"
opts <- parse_v2_args(require = character(0))
allow_unlocked <- isTRUE(opts$allow_unlocked)

cfg <- load_config("scz_ea_targeted")
scz <- load_config("schizophrenia")
assert_locked(list(scz_ea_targeted = cfg, schizophrenia = scz), allow_unlocked = allow_unlocked)
regions <- names(cfg$accepted_runs)

check_accepted <- function(module, cohort, re, want) {
    up <- require_accepted_upstream(module, cohort, re, allow_unaccepted = allow_unlocked)
    if (!identical(as.character(up$run_id), want)) {
        stop("config names ", want, " for ", module, " ", cohort, " ", re,
             " but the accepted run is ", up$run_id)
    }
    up
}
up_scz <- lapply(regions, function(re) check_accepted(MODULE, "AA", re, cfg$accepted_runs[[re]]))
up_cmb <- lapply(regions, function(re)
    check_accepted("05_cpg_meqtl_burden", "AA", re, cfg$meqtl_runs[[re]]))
names(up_scz) <- names(up_cmb) <- regions
for (re in regions) {
    for (p in c(file.path("01b_estimation_cells", "_m", "runs", cfg$estimation_cells[[re]]),
                file.path("01_vmr_catalog", "_m", "runs", cfg$catalog_runs[[re]]))) {
        if (!dir.exists(file.path(repo_root(), p))) stop("Missing upstream run ", p)
    }
}

scz_dir <- function(re) file.path(repo_root(), MODULE, "_m", "runs", cfg$accepted_runs[[re]])
targets <- rbindlist(lapply(regions, function(re) {
    pl <- fread(file.path(scz_dir(re), "results", "prioritized-loci.tsv"))[prioritized == TRUE]
    if (nrow(pl) > scz$prioritization$max_loci) stop(re, ": more prioritized loci than max_loci")
    tt <- fread(cmd = paste("zcat", shQuote(file.path(scz_dir(re), "results",
                                                     "risk-variant-cpg-tests.tsv.gz"))))
    tt <- tt[significant == TRUE & locus_id %in% pl$locus_id]
    tt[, .(region = re, locus_id, index_snp, risk_variant_id, ld_r2, cpg_id, vmr_id, chrom,
           aa_slope = slope, aa_slope_se = slope_se, aa_pval = pval_nominal, aa_qvalue = qvalue,
           aa_af = af, start_distance)]
}))
if (!nrow(targets)) stop("No nominated pairs at the prioritized loci")
setorder(targets, region, locus_id, aa_pval, risk_variant_id, cpg_id)
targets[, is_lead_pair := seq_len(.N) == 1L, by = .(region, locus_id)]

loci <- rbindlist(lapply(regions, function(re) {
    pl <- fread(file.path(scz_dir(re), "results", "prioritized-loci.tsv"))[prioritized == TRUE]
    pl[, .(region = re, rank, locus_id, chrom, start, end, index_snp_prioritized = index_snp,
           pgc3_index_pvalue)]
}))
loci <- merge(loci, targets[is_lead_pair == TRUE,
                            .(region, locus_id, lead_risk_variant_id = risk_variant_id,
                              lead_cpg_id = cpg_id, lead_aa_slope = aa_slope, lead_aa_pval = aa_pval)],
              by = c("region", "locus_id"), all.x = TRUE)
loci <- merge(loci, targets[, .(n_nominated_pairs = .N, n_cpgs = uniqueN(cpg_id),
                                n_variants = uniqueN(risk_variant_id)), by = .(region, locus_id)],
              by = c("region", "locus_id"), all.x = TRUE)
setorder(loci, region, rank)

if (!is.null(opts$run_id) && !allow_unlocked) stop("--run-id may only be given for a smoke run")
run <- new_run(
    module = "scz", cohort = COHORT, region = "crossregion",
    module_root = file.path(repo_root(), MODULE), run_id = opts$run_id,
    upstream = c(
        stats::setNames(lapply(regions, function(re) cfg$accepted_runs[[re]]), paste0("scz_aa_", regions)),
        stats::setNames(lapply(regions, function(re) cfg$meqtl_runs[[re]]), paste0("meqtl_aa_", regions)),
        stats::setNames(lapply(regions, function(re) cfg$estimation_cells[[re]]), paste0("estcell_ea_", regions)),
        stats::setNames(lapply(regions, function(re) cfg$catalog_runs[[re]]), paste0("vmrcat_pooled_", regions))),
    extra = list(
        smoke_run = if (allow_unlocked) "TRUE" else "FALSE",
        analysis = "ea_targeted_illustrative_loci",
        gating = "FALSE",
        config_scz_ea_targeted_sha256 = attr(cfg, "config_sha256"),
        config_schizophrenia_sha256 = attr(scz, "config_sha256"),
        n_loci = nrow(loci), n_pairs = nrow(targets)))
for (d in c("results", "work", "code/config")) dir.create(file.path(run$dir, d), recursive = TRUE, showWarnings = FALSE)
invisible(file.copy(file.path(repo_root(), "config",
                              c("scz_ea_targeted.yml", "schizophrenia.yml", "covariates.yml",
                                "meqtl_parameters.yml", "paths.yml")),
                    file.path(run$dir, "code", "config")))
write_atomic(targets, file.path(run$dir, "results", "targets.tsv"))
write_atomic(loci, file.path(run$dir, "results", "loci.tsv"))
print(loci[, .(region, rank, locus_id, index_snp_prioritized, n_nominated_pairs, lead_aa_pval)])
message("[09.19] opened ", run$run_id, ": ", nrow(loci), " loci, ", nrow(targets), " pairs")
cat(run$run_id, "\n", sep = "")
