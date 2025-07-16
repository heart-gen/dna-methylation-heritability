#### Calculate mean and SD of methylation values ####

suppressPackageStartupMessages({
    library('bsseq')
    library('HDF5Array')
    library(DelayedMatrixStats)
    library('data.table')
    library(here)
    library(dplyr)
    library(genio)
    library(plinkr)
})

args <- commandArgs(trailingOnly = TRUE)
chr  <- args[1]

## Function
filter_pheno <- function(BSobj, pheno_file_path) {
    pheno <- fread(pheno_file_path, header = TRUE)
    pheno_filtered <- pheno %>%
        filter(race == "AA", agedeath >= 17, region == "DLPFC")
    id    <- intersect(pheno_filtered$brnum, colData(BSobj)$brnum)
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

extract_fid_iid <- function(psam_file) {
    samples <- read_plink2_psam_file(psam_file)
    samples <- samples[, 1:2]
    return(samples)
}

write_meth_to_phen <- function(BSobj, M, samples, out_cpg) {
  
                                        # get CpG positions
    rownames(M) <- start(BSobj)
    
                                        # add sample IDs to methylation matrix
    sample_ids <- colData(BSobj)$brnum
    meth_df <- data.frame(FID = sample_ids, t(M), check.names = FALSE)

                                        # merge methylation data and sample ids by FID
    meth_merged <- meth_df %>%
        inner_join(samples, by = "FID") %>%
        arrange(match(FID, samples$FID)) %>%
        select(FID, IID, everything())
    
                                        # write methylation values to .phen file
    f_p <- file.path(out_cpg, "cpg_meth.phen")
    fwrite(meth_merged, file = f_p, sep = "\t", col.names = TRUE)

                                        # write CpG names
    f_p_names <- file.path(out_cpg, "cpg_pos.txt")
    write.table(colnames(meth_merged), file = f_p_names, row.names = FALSE,
                col.names = FALSE, quote = FALSE)
    
    return(meth_merged)
}

write_covar <- function(BSobj, pheno, id, meth_merged, out_covs) {
    out_cov  <- file.path(out_covs, "TOPMed_LIBD.AA.covar")
    out_qcov <- file.path(out_covs, "TOPMed_LIBD.AA.qcovar")
                                        # Filter data
    filtered_pheno <- pheno |>
        select(brnum, sex, primarydx, agedeath) |> filter(brnum %in% id)
                                        # Merge data
    meth_selected <- meth_merged |>
        inner_join(filtered_pheno, by= c("FID" = "brnum")) |>
        arrange(match(FID, meth_merged$FID)) |>
        tibble::column_to_rownames("FID")
                                        # Write file
    covar_merged <- meth_selected |> select(IID, sex, primarydx)
    covar_merged |>
        write.table(file=out_cov, sep="\t", row.names=TRUE,
                    col.names=FALSE, quote=FALSE)
    qcovar_merged <- meth_selected |> select(IID, agedeath)
    qcovar_merged |>
        write.table(file=out_qcov, sep="\t", row.names=TRUE,
                    col.names=FALSE, quote=FALSE)
  return(list(covar_merged=covar_merged, qcovar_merged=qcovar_merged))
}

## Main
                                        # load data
load(here("inputs/wgbs-data/dlpfc", paste0("dlpfc_chr", chr, 
    "BSobj_GenotypesDosage.rda")))
output_path <- here("heritability", "dlpfc", "_m")
subdirs <- c("covs", "cpg")

                                        # create output directories if they  
                                        # don't exist
for (subdir in subdirs) {
    subdir_path <- file.path(output_path, subdir, paste0("chr_", chr))
    if (!dir.exists(subdir_path)) {
        dir.create(subdir_path, recursive = TRUE)
    }
}

                                        # define output directories 
out_covs  <- file.path(output_path, "covs",  paste0("chr_", chr))
out_cpg   <- file.path(output_path, "cpg", paste0("chr_", chr))

                                        # keep only adult AA
pheno_file_path <- here("inputs/phenotypes/_m/phenotypes-AA.tsv")
filtered        <- filter_pheno(BSobj, pheno_file_path)

                                        # exclude low coverage sites
BSobj <- exclude_low_cov(filtered$BSobj)

                                        # calculate sd & mean of DNAm
stats <- DNAm_stats(BSobj, file.path(out_cpg, "stats.rda"))

                                        # read in FID, IID from sample file
psam_file <- here("inputs/genotypes/TOPMed_LIBD.AA.psam")
samples   <- extract_fid_iid(psam_file)

                                        # merge methylation values with FID and
                                        # IID and write to .phen file
meth_merged <- write_meth_to_phen(BSobj, stats$M, samples, out_cpg)

                                        # write covariate files
covars <- write_covar(BSobj, filtered$pheno, filtered$id, meth_merged, out_covs)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
