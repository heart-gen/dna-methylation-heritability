#### Is the caudate VMR/CpG universe advantage sample size, or coverage depth? ####
##
## Module 10 downsampled donors for meQTL *mapping* but inherited the VMR/CpG
## universe from full-N caudate. The universe itself is therefore still
## confounded, and it is the universe -- not the mapping -- that drives the
## common-universe rate comparison that failed the caudate gate.
##
## The universe is built by two steps in vmr-analysis/<region>/_h:
##
##   01.get_cpg_stats.R  exclude_low_cov():
##       keep CpGs with rowSums2(coverage >= 5) >= n_donors * 0.8
##
##   02c.write_top_cpg.R:
##       VMR seeds = CpGs with residual sd above the *within-region, per-chr
##       99th percentile* -- a RELATIVE cutoff, so seed count is mechanically
##       ~1% of whatever survives the coverage filter.
##
## Observed (all autosomes, full N):
##   caudate      24,435,609 CpGs pass QC -> 244,367 seeds   (N=153)
##   dlpfc        21,115,720                211,168          (N=111)
##   hippocampus  21,383,452                213,845          (N=116)
##
## Caudate keeps 15.7% more CpGs and gets 15.7% more seeds by construction.
## That tracks the VMR surplus (11,373 vs 9,976 / 9,801 = +14% / +16%) almost
## exactly, so step 1 is where the whole difference is created.
##
## Two candidate causes, which this script separates:
##   (a) SAMPLE SIZE. The 80%-of-donors threshold is a fixed fraction, but at
##       n=111 a CpG whose true adequate-coverage proportion sits near 0.8 is
##       sampled more noisily than at n=153.
##   (b) COVERAGE DEPTH. Caudate libraries may simply be sequenced deeper, so
##       more CpGs genuinely clear cov>=5 in >=80% of donors. That is a
##       sequencing property, not sample size and not biology.
##
## Design: re-run ONLY the coverage filter, on the same 30 N=111 caudate donor
## subsets module 10 used (seed 20260805), and compare against full-N caudate
## and against the two comparator regions. Also record per-donor coverage so
## (b) is measured directly, and recompute the residual-sd 99th percentile on
## the N-matched universe to show which absolute variability bar each region's
## "top 1%" actually corresponds to.
##
## Reading the output:
##   caudate N-matched ~= caudate full-N, still > comparators -> cause (b)
##   caudate N-matched falls to comparator level                -> cause (a)

suppressPackageStartupMessages({
    library('bsseq')
    library('HDF5Array')
    library('DelayedMatrixStats')
    library('data.table')
    library('matrixStats')
    library('here')
    library('dplyr')
})

args   <- commandArgs(trailingOnly = TRUE)
chr    <- args[1]
region <- args[2]

PROJECT   <- here()
OUTDIR    <- file.path(PROJECT, "meqtl-validation/13_vmr_universe_nmatched/_m")
REP_LISTS <- file.path(PROJECT, "meqtl-validation/04_cross_region_sharing/_m",
                       "caudate_downsample/sample_lists")
MIN_COV   <- 5      # matches exclude_low_cov()
MIN_FRAC  <- 0.8    # matches exclude_low_cov()
SD_QUANT  <- 0.99   # matches 02c.write_top_cpg.R
N_PC      <- 5      # matches 02b.res_var.R

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ---- helpers mirroring vmr-analysis/<region>/_h/01.get_cpg_stats.R ----------

filter_pheno <- function(BSobj, pheno_file_path, region) {
    pheno <- fread(pheno_file_path, header = TRUE)
    pheno_filtered <- pheno %>%
        filter(race == "AA", agedeath >= 17, region == !!region)
    id    <- intersect(pheno_filtered$brnum, colData(BSobj)$brnum)
    BSobj <- BSobj[, colData(BSobj)$brnum %in% id]
    return(list(BSobj = BSobj, pheno = pheno_filtered, id = id))
}

remove_ct_snps <- function(f_snp, BSobj) {
    snp <- fread(f_snp, header = FALSE, data.table = FALSE)[, 1]
    idx <- is.element(start(BSobj), snp)
    return(BSobj[!idx, ])
}

## The dlpfc and hippocampus BSobj archives have broken HDF5 backing:
##
##   M, Cov  point at <region>/_m/combined_hdf5/<region>_assays.h5, but the
##           file on disk is combined_hdf5/assays.h5 -- a rename, recoverable.
##   coef    points at an R session temp dump (auto<hash>.h5) that no longer
##           exists -- unrecoverable.
##
## Caudate is unaffected: all three of its assays live in one CpGassays.h5 that
## is still present. Any access, including the colData read inside filter_pheno,
## touches the seeds and fails.
##
## `coef` holds bsseq smoothing coefficients and is dropped rather than chased:
## this analysis needs only M and Cov, and dropping it makes hasBeenSmoothed()
## FALSE so getMeth() returns raw M/Cov everywhere. That is the point -- a
## between-region comparison must compute methylation the same way in all three
## regions, and smoothed values are only available for caudate.
drop_broken_coef <- function(BSobj) {
    if (!"coef" %in% assayNames(BSobj)) return(BSobj)
    a <- assay(BSobj, "coef", withDimnames = FALSE)
    s <- tryCatch(DelayedArray::seed(a), error = function(e) NULL)
    if (!is.null(s) && is(s, "HDF5ArraySeed") && !file.exists(s@filepath)) {
        message("dropping unrecoverable smoothing coefficients: ", s@filepath)
        keep <- setdiff(assayNames(BSobj), c("coef", "se.coef"))
        assays(BSobj, withDimnames = FALSE) <-
            assays(BSobj, withDimnames = FALSE)[keep]
    }
    BSobj
}

repoint_h5 <- function(BSobj, region) {
    a0 <- assay(BSobj, 1L, withDimnames = FALSE)
    if (!is(a0, "DelayedArray")) return(BSobj)
    stale <- DelayedArray::seed(a0)@filepath
    if (file.exists(stale)) return(BSobj)

    correct <- file.path("/projects/b1213/resources/libd_data/wgbs/new-data",
                         region, "_m", "combined_hdf5", "assays.h5")
    if (!file.exists(correct)) {
        stop("HDF5 backing file not found for ", region,
             "\n  seed points at: ", stale,
             "\n  tried instead:  ", correct)
    }
    message("repointing HDF5 backing: ", stale, " -> ", correct)

    for (nm in assayNames(BSobj)) {
        a <- assay(BSobj, nm, withDimnames = FALSE)
        a <- DelayedArray::modify_seeds(a, function(s) {
            if (is(s, "HDF5ArraySeed")) s@filepath <- correct
            s
        })
        assay(BSobj, nm, withDimnames = FALSE) <- a
    }
    BSobj
}

## ---- load ------------------------------------------------------------------

message("region=", region, " chr=", chr)
load(file.path(PROJECT, "inputs/wgbs-data", region, "_m",
               paste0(region, "_chr", chr, "_BSobj.rda")))

BSobj      <- repoint_h5(drop_broken_coef(BSobj), region)
pheno_file <- file.path(PROJECT, "inputs/phenotypes/_m/phenotypes-AA.tsv")
filtered   <- filter_pheno(BSobj, pheno_file, region)
BSobj      <- filtered$BSobj

if (!chr %in% c("X", "Y")) {
    f_snp <- paste0("/projects/b1213/resources/libd_data/wgbs/DEM2/snps_CT/chr", chr)
    BSobj <- remove_ct_snps(f_snp, BSobj)
}

donors  <- as.character(colData(BSobj)$brnum)
n_full  <- length(donors)
n_sites <- nrow(BSobj)
message("donors=", n_full, " sites_pre_cov_filter=", n_sites)

## ---- coverage matrix, thresholded once -------------------------------------
## `adequate` is CpG x donor logical; every downstream count is a column
## subset of it, so the expensive comparison happens exactly once.

cov      <- as.matrix(getCoverage(BSobj))
mean_cov <- colMeans2(cov)                 # taken before `cov` is released
adequate <- cov >= MIN_COV
rm(cov); gc()

## per-donor coverage, to measure cause (b) directly
donor_cov <- data.table(
    region      = region,
    chr         = chr,
    brnum       = donors,
    mean_cov    = mean_cov,
    frac_ge_min = colMeans2(adequate)
)
fwrite(donor_cov,
       file.path(OUTDIR, sprintf("donor_coverage.%s.chr%s.tsv", region, chr)),
       sep = "\t")

## ---- pass counts -----------------------------------------------------------

pass_count <- function(idx) {
    ## number of CpGs clearing the filter using only donors in `idx`
    n <- length(idx)
    sum(rowSums2(adequate, cols = idx) >= n * MIN_FRAC)
}

results <- list()

results[[length(results) + 1L]] <- data.table(
    region = region, chr = chr, design = "full_n", replicate = NA_integer_,
    n_donors = n_full, n_sites_pre = n_sites, n_sites_pass = pass_count(seq_len(n_full))
)

## N-matched arm: caudate only, reusing module 10's 30 replicate donor lists
if (region == "caudate") {
    rep_files <- sort(list.files(REP_LISTS, pattern = "^caudate_downsample_rep[0-9]+\\.txt$",
                                 full.names = TRUE))
    message("N-matched replicates found: ", length(rep_files))
    for (rf in rep_files) {
        rep_id  <- as.integer(sub(".*rep0*([0-9]+)\\.txt$", "\\1", rf))
        keep    <- readLines(rf)
        idx     <- which(donors %in% keep)
        if (length(idx) < 2) {
            warning("replicate ", rep_id, " matched ", length(idx), " donors; skipping")
            next
        }
        results[[length(results) + 1L]] <- data.table(
            region = region, chr = chr, design = "n_matched", replicate = rep_id,
            n_donors = length(idx), n_sites_pre = n_sites, n_sites_pass = pass_count(idx)
        )
    }
}

pass_dt <- rbindlist(results)
fwrite(pass_dt,
       file.path(OUTDIR, sprintf("coverage_pass.%s.chr%s.tsv", region, chr)),
       sep = "\t")

## ---- residual-sd cutoff, caudate only --------------------------------------
## 02c takes the top 1% *within* whatever survived the coverage filter, so the
## absolute variability bar is region-specific. The question this arm answers is
## narrow: does caudate's own bar move when its donors are N-matched?
##
## CAUDATE ONLY, deliberately. The bar is computed from smoothed methylation,
## and the dlpfc/hippocampus smoothing coefficients are unrecoverable (see
## drop_broken_coef). Recomputing their bars from raw M/Cov would compare a
## smoothed caudate value against unsmoothed comparators -- worse than not
## computing them. The cross-region bar comparison instead uses the values the
## pipeline itself already wrote, in
## vmr-analysis/<region>/_m/cpg/top1_cpg.tsv.

if (region != "caudate") {
    message("skipping sd-cutoff arm: smoothing coefficients unavailable for ", region)
    message("cross-region bars come from vmr-analysis/<region>/_m/cpg/top1_cpg.tsv")
    print(pass_dt)
    cat("\nReproducibility information:\n")
    print(Sys.time()); print(proc.time())
    options(width = 120)
    print(sessioninfo::session_info())
    quit(save = "no", status = 0)
}

meth <- as.matrix(getMeth(BSobj))
if (anyNA(meth)) {
    stop("unexpected NA in smoothed methylation for caudate chr", chr)
}
gc()

## PCA on a random CpG subsample. The PCs are a property of the donors, not of
## any particular CpG, so a subsample estimates them adequately -- and a full
## prcomp on a whole chromosome exhausts memory.
PCA_SITES <- 50000L
set.seed(20260805)

sd_cutoff_for <- function(idx) {
    n    <- length(idx)
    keep <- which(rowSums2(adequate, cols = idx) >= n * MIN_FRAC)
    if (length(keep) < 1000) return(NA_real_)

    pca_rows <- if (length(keep) > PCA_SITES) sample(keep, PCA_SITES) else keep
    pcs <- prcomp(t(meth[pca_rows, idx, drop = FALSE]),
                  center = TRUE, scale. = FALSE)$x[, seq_len(min(N_PC, n - 1)), drop = FALSE]
    design <- cbind(1, pcs)
    hat    <- diag(n) - design %*% solve(crossprod(design), t(design))

    ## residualise and accumulate row sds in chunks to bound peak memory
    sds   <- numeric(length(keep))
    chunk <- 100000L
    for (s in seq(1, length(keep), by = chunk)) {
        e   <- min(s + chunk - 1L, length(keep))
        blk <- meth[keep[s:e], idx, drop = FALSE]
        sds[s:e] <- rowSds(blk %*% hat)
        rm(blk)
    }
    gc()
    as.numeric(quantile(sds, probs = SD_QUANT, na.rm = TRUE))
}

sd_rows <- list(data.table(
    region = region, chr = chr, design = "full_n", replicate = NA_integer_,
    n_donors = n_full, sd_cutoff = sd_cutoff_for(seq_len(n_full))
))

if (region == "caudate") {
    ## 5 replicates is enough to bound the cutoff; the PCA is the costly part
    rep_files <- sort(list.files(REP_LISTS, pattern = "^caudate_downsample_rep[0-9]+\\.txt$",
                                 full.names = TRUE))[1:5]
    for (rf in rep_files) {
        rep_id <- as.integer(sub(".*rep0*([0-9]+)\\.txt$", "\\1", rf))
        idx    <- which(donors %in% readLines(rf))
        sd_rows[[length(sd_rows) + 1L]] <- data.table(
            region = region, chr = chr, design = "n_matched", replicate = rep_id,
            n_donors = length(idx), sd_cutoff = sd_cutoff_for(idx)
        )
    }
}

fwrite(rbindlist(sd_rows),
       file.path(OUTDIR, sprintf("sd_cutoff.%s.chr%s.tsv", region, chr)),
       sep = "\t")

print(pass_dt)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
