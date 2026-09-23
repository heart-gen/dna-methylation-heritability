#### 01_vmr_catalog / 00_prepare: per-chromosome CpG matrix and covariates ####
##
## Replaces vmr-analysis/{,all_individuals/}{caudate,dlpfc,hippocampus}/_h/
## 01.get_cpg_stats.R -- six near-identical copies -- with one parameterized
## script (AGENTS.md 5.3).
##
## Usage:
##   Rscript 00_prepare.R --cohort AA --region dlpfc --chrom 21 --run-id ID
##
## Repairs applied here:
##   V2  The region filter comes from --region. The legacy BA_only copies all
##       hard-coded `region == "caudate"`, so DLPFC and hippocampus discovered
##       VMRs on ~15 fewer donors than they were later modeled on.
##   V3  No here() calls to a stale `heritability/<region>` root; every path is
##       resolved from config/paths.yml.
##   V4  Sex chromosomes have no C->T mask and are not part of the primary
##       catalog. They are written under excluded/ with an explicit manifest.
##   V5  remove_ct_snps() operates on its argument, not on a global.
##   V8  read_psam() handles the headerless AA .psam.
##
## F14 Every donor and CpG exclusion decided here is now recorded, not only
##     messaged. 04_turnover.R assembles the per-chromosome parts into
##     qc/exclusions.tsv and qc/exclusion_accounting.tsv (AGENTS.md 7.1, 11).
##     The filters themselves are unchanged: the survivor set is asserted
##     against the ledger below, so recording cannot alter the catalog.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(V2_ROOT, "01_vmr_catalog", "_h", "exclusion_ledger.R"))

suppressPackageStartupMessages({
    library(bsseq)
    library(HDF5Array)
    library(DelayedMatrixStats)
    library(data.table)
    library(dplyr)
})

opts <- parse_v2_args(require = c("cohort", "region", "chrom", "run_id"))
cohort <- opts$cohort; region <- opts$region; chrom <- as.character(opts$chrom)

th   <- load_config("thresholds")
covs <- load_config("covariates")
assert_locked(list(thresholds = th, cohorts = load_config("cohorts")),
              allow_unlocked = opts$allow_unlocked)

arm <- cohort_def(cohort)
vmr_cfg <- th$vmr

## ---------------------------------------------------------------- functions

#' Donors for this cohort x region, from the phenotype table.
#'
#' The region filter is `region == region_arg`. Not a literal. This one line is
#' defect V2.
#'
#' Returns the kept donors plus the material an exclusion ledger needs: the
#' candidate set this region's catalog was drawn from, and one reason per
#' candidate that did not survive. `keep` is derived by the same expression as
#' before, from the same table, so the ledger observes the filter rather than
#' reimplementing it -- and the caller asserts the two agree.
select_donors <- function(pheno_file, region_arg, race_filter, min_age,
                          blacklist = NULL) {
    pheno <- fread(pheno_file, header = TRUE)
    for (col in c("brnum", "region", "race", "agedeath")) {
        if (!col %in% names(pheno)) {
            stop("Phenotype table is missing required column '", col, "': ",
                 pheno_file)
        }
    }
    keep <- pheno[agedeath >= min_age &
                  region == region_arg &
                  race %in% race_filter]
    if (!is.null(blacklist) && length(blacklist) > 0) {
        n_before <- nrow(keep)
        keep <- keep[!brnum %in% blacklist]
        message("[donors] blacklist removed ", n_before - nrow(keep), " donor(s)")
    }
    assert_no_dups(keep$brnum, "brnum in phenotype table")

    ## Candidate set: the phenotype rows for THIS brain region. Rows for other
    ## regions are not donors this catalog could have had, so they are reported
    ## as an aggregate rather than as excluded donors -- the same person is a
    ## candidate in the other region's run.
    in_region <- !is.na(pheno$region) & pheno$region == region_arg
    cand <- pheno[in_region]
    ## Precedence order fixes which single rule a multiply-failing donor is
    ## counted under, so the ledger's counts add up.
    reasons <- first_reason(list(
        donor_below_min_age_or_age_missing =
            is.na(cand$agedeath) | cand$agedeath < min_age,
        donor_race_not_in_cohort_arm =
            is.na(cand$race) | !(cand$race %in% race_filter),
        donor_on_sample_blacklist =
            cand$brnum %in% (blacklist %||% character())))

    list(keep = keep, candidates = cand, reasons = reasons,
         n_rows_other_region = sum(!in_region))
}

#' Drop CpGs overlapping C->T SNPs, which masquerade as unmethylated cytosines.
#'
#' V5: the legacy version read `filtered$BSobj` from the enclosing scope and
#' ignored its own BSobj argument.
remove_ct_snps <- function(BSobj, ct_file) {
    snp <- fread(ct_file, header = FALSE, data.table = FALSE)[, 1]
    idx <- is.element(start(BSobj), snp)
    message("[ct] removed ", sum(idx), " CpG(s) overlapping C->T SNPs")
    BSobj[!idx, ]
}

#' Keep CpGs covered at >= min_coverage in >= min_covered_fraction of donors.
exclude_low_cov <- function(BSobj, min_coverage, min_fraction) {
    cov <- getCoverage(BSobj)
    n <- ncol(BSobj)
    keep <- which(rowSums2(cov >= min_coverage) >= n * min_fraction)
    message("[cov] kept ", length(keep), " of ", nrow(BSobj), " CpGs at >=",
            min_coverage, "x in >=", min_fraction * 100, "% of ", n, " donors")
    BSobj[keep, ]
}

## ------------------------------------------------------------------- main

module_root <- file.path(V2_ROOT, "01_vmr_catalog")
run_dir <- file.path(module_root, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) {
    stop("Run directory not found: ", run_dir,
         "\n  Create it with step_1.sh, which calls new_run() once per ",
         "cohort x region before fanning out over chromosomes.")
}

is_primary <- is_primary_chrom(chrom)
out_base <- if (is_primary) run_dir else file.path(run_dir, "excluded")
out_cpg  <- file.path(out_base, "cpg", paste0("chr_", chrom))
out_covs <- file.path(out_base, "covs", paste0("chr_", chrom))
dir.create(out_cpg, recursive = TRUE, showWarnings = FALSE)
dir.create(out_covs, recursive = TRUE, showWarnings = FALSE)

if (!is_primary) {
    message("[chrom] chr", chrom, " is a sex chromosome: writing to excluded/. ",
            "It has no C->T mask and is not part of the primary catalog (V4).")
}

## Load WGBS
BSobj <- load_bsobj(region, chrom)
bsobj_source <- attr(BSobj, "v2_source")
bsobj_sha256 <- attr(BSobj, "v2_source_sha256")

## Donors
blacklist <- sample_blacklist(region)
sel <- select_donors(
    pheno_file   = arm$phenotype_table,
    region_arg   = region,
    race_filter  = unlist(arm$race_filter),
    min_age      = vmr_cfg$min_age,
    blacklist    = blacklist
)
pheno <- sel$keep

samples <- read_psam(arm$psam)

## The analysis set is the intersection of three sources. Taking it explicitly,
## in one place, is what lets assert_expected_n() mean anything.
bs_ids <- as.character(colData(BSobj)$brnum)
assert_no_dups(bs_ids, "brnum in BSobj colData")
analysis_ids <- Reduce(intersect, list(pheno$brnum, bs_ids, samples$FID))
## Deterministic order: phenotype-table order, filtered.
analysis_ids <- pheno$brnum[pheno$brnum %in% analysis_ids]

if (length(analysis_ids) == 0) {
    stop("No donors survive the intersection of phenotype table, BSobj, and ",
         "psam for ", cohort, "/", region)
}
message("[donors] phenotype ", nrow(pheno), " | BSobj ", length(bs_ids),
        " | psam ", nrow(samples), " -> analysis set ", length(analysis_ids))
assert_expected_n(length(analysis_ids), cohort, region)

## ------------------------------------------------- donor exclusion ledger (F14)
##
## The two input-availability rules are applied to the candidate set here, after
## the same two `intersect` terms that produced analysis_ids above, so the ledger
## reads off the filter instead of restating it.
cand <- sel$candidates
donor_reasons <- sel$reasons
avail <- list(donor_absent_from_wgbs_bsobj = !(cand$brnum %in% bs_ids),
              donor_absent_from_genotype_psam = !(cand$brnum %in% samples$FID))
for (nm in names(avail)) {
    hit <- is.na(donor_reasons) & avail[[nm]]
    donor_reasons[hit] <- nm
}

## The property that makes this a record and not a second filter: the donors the
## ledger leaves unexcluded are exactly the donors the pipeline kept. If this
## ever fails, the ledger and the pipeline disagree and the run must stop --
## a wrong exclusion table is worse than none.
ledger_survivors <- cand$brnum[is.na(donor_reasons)]
if (!setequal(ledger_survivors, analysis_ids) ||
    length(ledger_survivors) != length(analysis_ids)) {
    stop("Donor exclusion ledger does not reproduce the analysis set: ",
         length(ledger_survivors), " ledger survivors against ",
         length(analysis_ids), " analysis donors.\n  Ledger-only: ",
         paste(head(setdiff(ledger_survivors, analysis_ids), 5), collapse = ", "),
         " | analysis-only: ",
         paste(head(setdiff(analysis_ids, ledger_survivors), 5), collapse = ", "))
}

donor_ledger <- rbindlist(list(
    excl_rows(stage = "00_prepare", unit_type = "donor",
              exclusion_reason = donor_reasons[!is.na(donor_reasons)],
              unit_id = cand$brnum[!is.na(donor_reasons)],
              chrom = chrom),
    ## Not listable as donors of this catalog: the same people are candidates in
    ## the run for the region their row names.
    excl_rows(stage = "00_prepare", unit_type = "phenotype_row",
              exclusion_reason = "phenotype_row_for_another_brain_region",
              unit_id = NA_character_, chrom = chrom,
              n_units = sel$n_rows_other_region, itemized = FALSE)),
    use.names = TRUE)
if (!is.null(blacklist) && length(blacklist) == 0) {
    message("[donors] no sample blacklist for ", region,
            " (retired in v2); recorded as zero exclusions")
}
write_atomic(donor_ledger, file.path(out_cpg, "exclusions_donors.tsv"))

BSobj <- BSobj[, match(analysis_ids, bs_ids)]
stopifnot(identical(as.character(colData(BSobj)$brnum), analysis_ids))

## CpG filtering. Counts are taken around the existing calls rather than inside
## them, so the filters are untouched and the arithmetic is observable:
## n_cpgs_input = n_cpgs + ct removed + low-coverage removed.
n_cpgs_input <- nrow(BSobj)
if (has_ct_mask(chrom)) {
    BSobj <- remove_ct_snps(BSobj, resolve_path("ct_snp_template", chrom = chrom,
                                                check = TRUE))
} else {
    message("[ct] no C->T mask for chr", chrom, "; CpGs are UNMASKED")
}
n_after_ct <- nrow(BSobj)
BSobj <- exclude_low_cov(BSobj, vmr_cfg$min_coverage, vmr_cfg$min_covered_fraction)

n_ct_removed <- n_cpgs_input - n_after_ct
n_lowcov_removed <- n_after_ct - nrow(BSobj)

if (nrow(BSobj) == 0) {
    stop("No CpGs survive QC on chr", chrom, " for ", cohort, "/", region)
}

## CpG exclusions are aggregate rows: itemizing millions of positions per
## chromosome would dwarf the catalog it documents, and the counts are what a
## denominator question needs. A CpG unmasked because no C->T list exists for
## this chromosome is recorded too -- it is the sex-chromosome hazard V4 names,
## and a zero would otherwise be indistinguishable from "we did not check".
cpg_ledger <- rbindlist(list(
    excl_rows("00_prepare", "cpg", "cpg_overlaps_ct_snp", NA_character_, chrom,
              n_ct_removed, itemized = FALSE),
    excl_rows("00_prepare", "cpg",
              paste0("cpg_below_min_coverage_", vmr_cfg$min_coverage, "x_in_",
                     vmr_cfg$min_covered_fraction * 100, "pct_of_donors"),
              NA_character_, chrom, n_lowcov_removed, itemized = FALSE)),
    use.names = TRUE)
if (!has_ct_mask(chrom)) {
    cpg_ledger <- rbindlist(list(cpg_ledger, excl_rows(
        "00_prepare", "cpg", "no_ct_snp_mask_available_cpgs_unmasked",
        NA_character_, chrom, 0L, itemized = FALSE)), use.names = TRUE)
}
write_atomic(cpg_ledger, file.path(out_cpg, "exclusions_cpgs.tsv"))

## Methylation matrix and per-CpG summaries
M <- as.matrix(getMeth(BSobj))
rownames(M) <- start(BSobj)
colnames(M) <- analysis_ids
sds   <- rowSds(M)
means <- rowMeans2(M)
save(sds, means, BSobj, file = file.path(out_cpg, "stats.rda"))

## Write the CpG matrix in donor order, with FID/IID from the psam.
aligned <- align_by_id(t(M), samples, id_x = analysis_ids, id_y = samples$FID,
                       ids = analysis_ids)
meth_out <- data.table(FID = analysis_ids,
                       IID = aligned$y$IID,
                       aligned$x)
setnames(meth_out, c("FID", "IID", as.character(start(BSobj))))
write_atomic(meth_out, file.path(out_cpg, "cpg_meth.phen"))
write_atomic(names(meth_out), file.path(out_cpg, "cpg_pos.txt"))

## Binary twin of cpg_meth.phen for 01_analyze.R.
##
## The .phen file is donors x CpGs, so a large chromosome is ~120 rows by 1.8M
## columns. fread() segfaults nondeterministically above roughly 1M columns
## (chr1/2/5/6/7/10 died on the first full-scale DLPFC run while chr3 read
## fine), and the smoke run never saw it because chr22 has only 492k columns.
## The .phen file is kept for downstream text consumers; analysis reads this.
meth_rds <- file.path(out_cpg, "cpg_meth.rds")
tmp_rds <- paste0(meth_rds, ".tmp")
saveRDS(list(FID = analysis_ids,
             IID = as.character(aligned$y$IID),
             cpg = as.character(start(BSobj)),
             meth = as.matrix(aligned$x)),
        tmp_rds, compress = FALSE)
if (!file.rename(tmp_rds, meth_rds)) {
    stop("Could not finalize ", meth_rds)
}

## Covariates, in the same donor order as the methylation matrix.
covar_src <- pheno[match(analysis_ids, brnum)]
stopifnot(identical(as.character(covar_src$brnum), analysis_ids))
write_atomic(
    data.table(FID = analysis_ids, IID = aligned$y$IID,
               sex = covar_src$sex, primarydx = covar_src$primarydx),
    file.path(out_covs, paste0(arm$covar_prefix, ".covar")), col.names = FALSE)
write_atomic(
    data.table(FID = analysis_ids, IID = aligned$y$IID,
               agedeath = covar_src$agedeath),
    file.path(out_covs, paste0(arm$covar_prefix, ".qcovar")), col.names = FALSE)

## Donor manifest for this chromosome. The combine step checks these agree.
write_atomic(
    data.table(FID = analysis_ids, IID = aligned$y$IID,
               race = covar_src$race, region = covar_src$region,
               agedeath = covar_src$agedeath, sex = covar_src$sex,
               order_index = seq_along(analysis_ids)),
    file.path(out_cpg, "donors.tsv"))

## The n_donor_* and n_cpg_* fields are the denominators 04_turnover.R balances
## the exclusion ledger against; they are additions, so consumers that select
## fields by name are unaffected.
write_atomic(
    data.table(
        field = c("cohort", "region", "chrom", "is_primary_chrom",
                  "n_donors", "n_cpgs", "donor_checksum", "ct_masked",
                  "blacklist_n", "bsobj_source", "bsobj_sha256",
                  "n_donor_candidate_rows", "n_donor_unique_candidates",
                  "n_donors_excluded", "n_phenotype_rows_other_region",
                  "n_cpgs_input", "n_cpgs_excluded_ct_snp",
                  "n_cpgs_excluded_low_coverage"),
        value = c(cohort, region, chrom, is_primary,
                  length(analysis_ids), nrow(BSobj), donor_checksum(analysis_ids),
                  has_ct_mask(chrom), length(blacklist %||% character()),
                  bsobj_source, bsobj_sha256,
                  nrow(cand), length(unique(cand$brnum)),
                  sum(!is.na(donor_reasons)), sel$n_rows_other_region,
                  n_cpgs_input, n_ct_removed, n_lowcov_removed)),
    file.path(out_cpg, "prepare_summary.tsv"))

message("[done] chr", chrom, ": ", nrow(BSobj), " CpGs x ",
        length(analysis_ids), " donors")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
