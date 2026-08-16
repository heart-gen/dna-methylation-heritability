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

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

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
    keep
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
bsobj_file <- resolve_path("wgbs_bsobj_template", region = region, chrom = chrom,
                           check = TRUE)
message("[load] ", bsobj_file)
load(bsobj_file)  # provides BSobj

## Donors
blacklist <- sample_blacklist(region)
pheno <- select_donors(
    pheno_file   = arm$phenotype_table,
    region_arg   = region,
    race_filter  = unlist(arm$race_filter),
    min_age      = vmr_cfg$min_age,
    blacklist    = blacklist
)

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

BSobj <- BSobj[, match(analysis_ids, bs_ids)]
stopifnot(identical(as.character(colData(BSobj)$brnum), analysis_ids))

## CpG filtering
if (has_ct_mask(chrom)) {
    BSobj <- remove_ct_snps(BSobj, resolve_path("ct_snp_template", chrom = chrom,
                                                check = TRUE))
} else {
    message("[ct] no C->T mask for chr", chrom, "; CpGs are UNMASKED")
}
BSobj <- exclude_low_cov(BSobj, vmr_cfg$min_coverage, vmr_cfg$min_covered_fraction)

if (nrow(BSobj) == 0) {
    stop("No CpGs survive QC on chr", chrom, " for ", cohort, "/", region)
}

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

write_atomic(
    data.table(
        field = c("cohort", "region", "chrom", "is_primary_chrom",
                  "n_donors", "n_cpgs", "donor_checksum", "ct_masked",
                  "blacklist_n", "bsobj_sha256"),
        value = c(cohort, region, chrom, is_primary,
                  length(analysis_ids), nrow(BSobj), donor_checksum(analysis_ids),
                  has_ct_mask(chrom), length(blacklist %||% character()),
                  file_sha256(bsobj_file))),
    file.path(out_cpg, "prepare_summary.tsv"))

message("[done] chr", chrom, ": ", nrow(BSobj), " CpGs x ",
        length(analysis_ids), " donors")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
