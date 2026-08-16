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
## Three candidate causes, which this script separates:
##   (a) SAMPLE SIZE. The 80%-of-donors threshold is a fixed fraction, but at
##       n=111 a CpG whose true adequate-coverage proportion sits near 0.8 is
##       sampled more noisily than at n=153.
##   (b) COVERAGE DEPTH. Caudate libraries sequenced deeper, so more CpGs
##       clear cov>=5 at all.
##   (c) COVERAGE UNIFORMITY. Caudate libraries more consistent between donors.
##       The filter asks whether >=80% of donors clear cov>=5 at a CpG, so it
##       is broken by the low tail of the donor distribution, not by the mean:
##       a shallower but tighter region keeps more CpGs. (b) and (c) are
##       distinct and can point in opposite directions.
##
## Design: re-run ONLY the coverage filter, on the same 30 N=111 caudate donor
## subsets module 10 used (seed 20260805), and compare against full-N caudate
## and against the two comparator regions. Record per-donor mean coverage AND
## the per-donor fraction of adequately covered sites, so (b) and (c) are
## measured separately. Recompute the residual-sd 99th percentile per region
## and for N-matched caudate, to show which absolute variability bar each
## region's "top 1%" actually corresponds to.
##
## Reading the output:
##   caudate N-matched falls to comparator level                -> (a)
##   caudate N-matched ~= full-N, and caudate deeper            -> (b)
##   caudate N-matched ~= full-N, caudate shallower but tighter -> (c)

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

## ---- input selection -------------------------------------------------------
##
## The dlpfc and hippocampus `<region>_chr<N>_BSobj.rda` archives in this repo
## have broken HDF5 backing and cannot be loaded:
##
##   M, Cov  point at combined_hdf5/<region>_assays.h5; the file on disk is
##           combined_hdf5/assays.h5 -- a rename.
##   coef    points at an R session temp dump (auto<hash>.h5) that is gone.
##
## Alexis re-exported both regions with saveHDF5SummarizedExperiment(), which
## writes a self-contained directory (assays.h5 + se.rds) with all three assays
## resolving and hasBeenSmoothed() TRUE. All 22 autosomes are present for both.
## Prefer those; they are the repaired copy of the same objects (chr22 dlpfc is
## 598,862 x 176 either way).
##
## Caudate is unaffected -- its three assays live in one CpGassays.h5 that is
## still present -- and has no re-export, so it loads from the .rda as before.
STAFF_WGBS <- file.path(
    "/projects/b1213/users/alexis/projects/dna-methylation-heritability",
    "inputs/wgbs-data"
)

load_bsobj <- function(region, chr) {
    se_dir <- file.path(STAFF_WGBS, region, "_m",
                        paste0(region, "_chr", chr, "_BSobj"))
    if (dir.exists(se_dir) && file.exists(file.path(se_dir, "se.rds"))) {
        message("loading repaired HDF5SummarizedExperiment: ", se_dir)
        obj <- HDF5Array::loadHDF5SummarizedExperiment(se_dir)
        attr(obj, "source") <- se_dir
        return(obj)
    }

    rda <- file.path(PROJECT, "inputs/wgbs-data", region, "_m",
                     paste0(region, "_chr", chr, "_BSobj.rda"))
    message("loading local archive: ", rda)
    env <- new.env(parent = emptyenv())
    load(rda, envir = env)
    obj <- get("BSobj", envir = env)

    a0 <- assay(obj, 1L, withDimnames = FALSE)
    if (is(a0, "DelayedArray")) {
        fp <- DelayedArray::seed(a0)@filepath
        if (!file.exists(fp)) {
            stop("HDF5 backing missing for ", region, " chr", chr, ": ", fp,
                 "\nNo repaired export found at ", se_dir)
        }
    }
    attr(obj, "source") <- rda
    obj
}

## ---- load ------------------------------------------------------------------

message("region=", region, " chr=", chr)
BSobj <- load_bsobj(region, chr)

if (!isTRUE(hasBeenSmoothed(BSobj))) {
    stop("BSseq object for ", region, " chr", chr, " is not smoothed; ",
         "the sd-cutoff arm requires smoothed methylation in every region")
}
pheno_file <- file.path(PROJECT, "inputs/phenotypes/_m/phenotypes-AA.tsv")
filtered   <- filter_pheno(BSobj, pheno_file, region)
BSobj      <- filtered$BSobj

if (!chr %in% c("X", "Y")) {
    f_snp <- paste0("/projects/b1213/resources/libd_data/wgbs/DEM2/snps_CT/chr", chr)
    BSobj <- remove_ct_snps(f_snp, BSobj)
}

## Drop donors without genotype PCs, as vmr-analysis/<region>/_h/02.pca.R does
## (get_snp_pcs removes any sample with NA in snpPC1-10, and res_snp_pcs then
## subsets the methylation matrix to the survivors). Exactly one donor per
## region is affected, which is the whole of the 154/112/117 vs 153/111/116
## discrepancy against module 10's design summary. It matters beyond
## bookkeeping: the coverage filter asks whether >=80% of *these* donors clear
## cov>=5, so including an ungenotyped donor would shift the universe.
this_region <- region
geno_ok <- fread(pheno_file, header = TRUE)[
    race == "AA" & agedeath >= 17 & region == this_region,
    c("brnum", paste0("snpPC", 1:10)), with = FALSE
]
geno_ok <- geno_ok[complete.cases(geno_ok), brnum]

all_donors <- as.character(colData(BSobj)$brnum)
dropped    <- setdiff(all_donors, geno_ok)
if (length(dropped) > 0) {
    message("dropping ", length(dropped), " donor(s) without genotype PCs: ",
            paste(dropped, collapse = ", "))
    BSobj <- BSobj[, all_donors %in% geno_ok]
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

## ---- residual-sd cutoff -----------------------------------------------------
## 02c takes the top 1% *within* whatever survived the coverage filter, so the
## absolute variability bar is region-specific. Two things are measured here:
## how the bar differs between regions, and whether caudate's own bar moves
## when its donors are N-matched.
##
## Computed for all three regions on smoothed methylation. That is possible
## only because the dlpfc/hippocampus objects were re-exported with intact
## smoothing coefficients (see load_bsobj) -- the broken .rda archives would
## have forced a smoothed-vs-raw comparison, which would not be interpretable.
##
## The residualisation reproduces the pipeline's two-stage adjustment:
##   02.pca.R   regress smoothed methylation on snpPC1-3 (genotype ancestry
##              PCs from the phenotype table), then PCA the residuals
##   02b.res_var.R  regress on the first 5 of those methylation PCs, take
##              rowSds of what is left
## Skipping the snpPC stage changes the PCs and reverses the between-region
## ordering of the bar, so it is not optional for a cross-region comparison.

meth <- as.matrix(getMeth(BSobj))
if (anyNA(meth)) {
    stop("unexpected NA in smoothed methylation for ", region, " chr", chr)
}
gc()

## snpPC1-3 for the donors in this object, in column order
SNP_PCS <- 3L
snp_pc_tab <- fread(pheno_file, header = TRUE)[
    race == "AA" & agedeath >= 17 & region == this_region,
    c("brnum", paste0("snpPC", seq_len(SNP_PCS))), with = FALSE
]
snp_pc_tab <- snp_pc_tab[match(donors, brnum)]
if (anyNA(snp_pc_tab)) {
    stop("missing snpPC1-", SNP_PCS, " for ", sum(!complete.cases(snp_pc_tab)),
         " donors in ", region, " chr", chr)
}
snp_pc_mat <- as.matrix(snp_pc_tab[, -1L])

## PCA on a random CpG subsample. The PCs are a property of the donors, not of
## any particular CpG, so a subsample estimates them adequately -- and a full
## prcomp on a whole chromosome exhausts memory.
PCA_SITES <- 50000L
set.seed(20260805)

residualise_on <- function(m, design) {
    ## returns m with the column space of `design` projected out (donors are
    ## columns of m, so the projection acts from the right)
    hat <- diag(ncol(m)) - design %*% solve(crossprod(design), t(design))
    m %*% hat
}

sd_cutoff_for <- function(idx) {
    n    <- length(idx)
    keep <- which(rowSums2(adequate, cols = idx) >= n * MIN_FRAC)
    if (length(keep) < 1000) return(NA_real_)

    ## stage 1: remove genotype ancestry PCs, as 02.pca.R does
    snp_design <- cbind(1, snp_pc_mat[idx, , drop = FALSE])

    ## stage 2: PCA those residuals, then remove the leading N_PC of them
    pca_rows <- if (length(keep) > PCA_SITES) sample(keep, PCA_SITES) else keep
    pca_in   <- residualise_on(meth[pca_rows, idx, drop = FALSE], snp_design)
    pcs <- prcomp(t(pca_in), center = TRUE, scale. = FALSE)$x[
        , seq_len(min(N_PC, n - SNP_PCS - 2L)), drop = FALSE]
    rm(pca_in)

    design <- cbind(snp_design, pcs)
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
