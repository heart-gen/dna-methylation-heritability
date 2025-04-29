#### Calculate mean and SD of methylation values ####

suppressPackageStartupMessages({
    library('bsseq')
    library('HDF5Array')
    library(DelayedMatrixStats)
    library('data.table')
    library(scales)
    library(here)
    library(dplyr)
    library(genio)
    library(plinkr)
})

args <- commandArgs(trailingOnly = TRUE)
chr  <- args[1]

## Function
change_file_path <- function(BSobj, raw_assays) {
    var <- c("M", "Cov", "coef")
    for (assay in var) {
        BSobj@assays@data@listData[[assay]]@seed@seed@filepath <- raw_assays
    }
    return(BSobj)
}

filter_pheno <- function(BSobj, pheno_file_path) {
    pheno <- read.csv(pheno_file_path, header = TRUE)
    pheno_filtered <- pheno %>%
        filter(Race == "AA", Age >= 17, Region == "DLPFC")
    id    <- intersect(pheno_filtered$BrNum, colData(BSobj)$brnum)
    BSobj <- BSobj[, colData(BSobj)$brnum %in% id]
    return(list(BSobj = BSobj, pheno = pheno_filtered, id = id))
}

exclude_low_cov <- function(BSobj) {
    cov   <- getCoverage(BSobj)
    n     <- length(colData(BSobj)$brnum)
    keep  <- which(rowSums2(cov >= 5) >= n * 0.8)
    BSobj <- BSobj[keep, ]
    return(BSobj)
}

DNAm_stats <- function(BSobj, out_stats) {
    M     <- as.matrix(getMeth(BSobj))
    sds   <- rowSds(M)
    means <- rowMeans2(M)
    save(sds, means, BSobj, file=out_stats)
    return(list(M = M, sds = sds, means = means))
}

get_vmr <- function(v, out_vmr) {
    sdCut  <- quantile(v[, 3], prob = 0.99, na.rm = TRUE)
    vmrs   <- c()
    v2     <- v[v$chr==chr, ]
    v2     <- v2[order(v2$start), ]
    isHigh <- rep(0, nrow(v2))
    isHigh[v2$sd > sdCut] <- 1
    vmrs0  <- bsseq:::regionFinder3(isHigh, as.character(v2$chr), 
                                   v2$start, maxGap = 1000)$up
    vmr    <- vmrs0[vmrs0$n > 5,1:3]
    write.table(vmr, file=out_vmr, col.names=F,
                row.names=F, sep="\t", quote=F)
    return(vmr)
}

extract_fid_iid <- function(psam_file) {
    samples <- read_plink2_psam_file(psam_file)
    samples <- samples[, 1:2]
    return(samples)
}

write_meth_to_phen <- function(BSobj, M, samples, out_phen) {

                                        # add sample IDs to methylation matrix
    sample_ids <- colData(BSobj)$brnum
    meth_df <- data.frame(FID = sample_ids, t(M))

                                        # merge methylation data and sample ids by FID
    meth_merged <- meth_df %>%
        inner_join(samples, by = "FID") %>%
        arrange(match(FID, samples$FID)) %>%
        select(FID, IID, everything())
    colnames(meth_merged)[1:3] <- c("fam", "id", "pheno")

                                        # write methylation values to .phen file
    write_phen(file=out_phen, meth_merged)
    return(meth_merged)
}

write_covar <- function(BSobj, pheno, id, meth_merged, out_covs) {
    out_cov  <- file.path(out_covs, "TOPMed_LIBD.AA.covar")
    out_qcov <- file.path(out_covs, "TOPMed_LIBD.AA.qcovar")
                                        # Filter data
    filtered_pheno <- pheno |>
        select(BrNum, Sex, Dx, Age) |> filter(BrNum %in% id)
                                        # Merge data
    meth_selected <- select(meth_merged, fam, id) |>
        inner_join(filtered_pheno, by= c("fam" = "BrNum")) |>
        arrange(match(fam, meth_merged$fam)) |>
        tibble::column_to_rownames("fam")
                                        # Write file
    covar_merged <- meth_selected |> select(id, Sex, Dx)
    covar_merged |>
        write.table(file=out_cov, sep="\t", row.names=TRUE,
                    col.names=FALSE, quote=FALSE)
    qcovar_merged <- meth_selected |> select(id, Age)
    qcovar_merged |>
        write.table(file=out_qcov, sep="\t", row.names=TRUE,
                    col.names=FALSE, quote=FALSE)
  return(list(covar_merged=covar_merged, qcovar_merged=qcovar_merged))
}

## Main
                                        # load data
load(here("inputs/wgbs-data/dlpfc", paste0("dlpfc_chr", chr, 
    "BSobj_GenotypesDosage.rda")))
output_path <- here("heritability", "gcta", "dlpfc", "_m")
subdirs <- c("vmr", "covs", "cpg")

                                        # create output directories if they  
                                        # don't exist
for (subdir in subdirs) {
    subdir_path <- file.path(output_path, subdir, paste0("chr_", chr))
    if (!dir.exists(subdir_path)) {
        dir.create(subdir_path, recursive = TRUE)
    }
}

                                        # define output directories 
out_vmr   <- file.path(output_path, "vmr",   paste0("chr_", chr))
out_covs  <- file.path(output_path, "covs",  paste0("chr_", chr))
out_cpg   <- file.path(output_path, "cpg", paste0("chr_", chr))

                                        # change file path for raw data
raw_assays  <- here("inputs/wgbs-data/dlpfc/raw/CpGassays.h5")
BSobj       <- change_file_path(BSobj, raw_assays)

                                        # keep only adult AA
pheno_file_path <- here("inputs/phenotypes/merged/_m/merged_phenotypes.csv")
filtered        <- filter_pheno(BSobj, pheno_file_path)

                                        # exclude low coverage sites
BSobj <- exclude_low_cov(filtered$BSobj)

                                        # calculate sd & mean of DNAm
stats <- DNAm_stats(BSobj, file.path(out_cpg, "stats.rda"))

                                        # extract VMRs
v   <- data.frame(chr = chr, start = start(BSobj), sd = stats$sds)
vmr <- get_vmr(v, file.path(out_vmr, "vmr.bed"))

                                        # read in FID, IID from sample file
psam_file <- here("inputs/genotypes/TOPMed_LIBD.AA.psam")
samples   <- extract_fid_iid(psam_file)

                                        # merge methylation values with FID and
                                        # IID and write to .phen file
meth_merged <- write_meth_to_phen(BSobj, stats$M, samples, file.path(out_cpg, "cpg_meth.phen"))

                                        # write covariate files
covars <- write_covar(BSobj, filtered$pheno, filtered$id, meth_merged, out_covs)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
