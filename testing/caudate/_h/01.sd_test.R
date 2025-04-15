#### Calculate mean and SD of methylation values ####


library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')
library(scales)
library(here)
library(dplyr)
library(genio)
library(plinkr)

args <- commandArgs(trailingOnly = TRUE)
chr  <- args[1]

##here::i_am("testing/caudate/_h/01.sd_test.R")
##setwd(here("testing/caudate/_h"))
## Function
change_file_path <- function(BSobj,raw_assays) {
    var <- c("M", "Cov", "coef")
    for (assay in var) {
        BSobj@assays@data@listData[[assay]]@seed@seed@filepath <- raw_assays
    }
    return(BSobj)
}

filter_pheno <- function(BSobj,pheno_file_path) {
    pheno <- read.csv(pheno_file_path,header=T)
    id <- intersect(
        pheno$BrNum[
                  pheno$Race == "AA" &
                  pheno$Age >= 17 &
                  pheno$Region == "Caudate"
              ],
        colData(BSobj)$brnum
    )
    BSobj <- BSobj[,is.element(colData(BSobj)$brnum,id)]
    return(BSobj)
}

exclude_low_cov <- function(BSobj) {
    cov=getCoverage(BSobj)
    n     <- length(colData(BSobj)$brnum)
    keep  <- which(rowSums2(cov >= 5) >= n * 0.8)
    BSobj <- BSobj[keep,]
    return(BSobj)
}

DNAm_stats <- function(BSobj,output_stats) {
    M     <- as.matrix(getMeth(BSobj))
    sds   <- rowSds(M)
    means <- rowMeans2(M)
    save(sds,means,BSobj,file=output_stats)
    return(list(M = M, sds = sds, means = means))
}

get_vmr <- function(v) {
    sdCut <- quantile(v[,3], prob = 0.99, na.rm = TRUE)
    vmrs  <- c()
    v2    <- v[v$chr==chr,]
    v2    <- v2[order(v2$start),]
    isHigh <- rep(0, nrow(v2))
    isHigh[v2$sd > sdCut] <- 1
    vmrs0 <- bsseq:::regionFinder3(isHigh, as.character(v2$chr), v2$start, maxGap = 1000)$up
    vmr   <- vmrs0[vmrs0$n > 5,1:3]
    write.table(vmr,file=file.path(output_path,"vmr_test.bed"),col.names=F,
                row.names=F,sep="\t",quote=F)
    return(vmr)
}

extract_fid_iid <- function(psam_file) {
    samples <- read_plink2_psam_file(psam_file)
    samples <- samples[,1:2]
    return(samples)
}

write_meth_to_phen <- function(BSobj,M,samples,output) {

                                        # add sample IDs to methylation matrix
    sample_ids <- colData(BSobj)$brnum #use this or id?
    meth_transpose <- t(M)
    meth_df <- data.frame(FID = sample_ids, meth_transpose)

                                        # merge methylation data and sample ids by FID
    meth_merged <- meth_df %>%
        inner_join(samples, by = "FID") %>%
        arrange(match(FID, samples$FID)) %>%
        select(FID, IID, everything())
    colnames(meth_merged)[1:3] <- c("fam", "id", "pheno")

                                        # write methylation values to .phen file
    write_phen(file=file.path(output_path,"cpg_meth.phen"), meth_merged)
    return(meth_merged)
}
write_covar <- function(BSobj, pheno_file_path, meth_merged, output_path) {
    out_cov  <- file.path(output_path, "TOPMed_LIBD.AA.covar")
    out_qcov <- file.path(output_path, "TOPMed_LIBD.AA.qcovar")
                                        # Filter data
    pheno   <- filter_pheno(BSobj,pheno_file_path) |>
        select(BrNum, Sex, Dx, Age) |> filter(BrNum %in% id)
                                        # Merge data
    meth_selected <- select(meth_merged, fam, id) |>
        inner_join(pheno, by="fam"="BrNum") |>
        arrange(match(fam, meth_merged$fam)) |>
        tibble::column_to_rownames("BrNum")
                                        # Write file
    covar_merged <- meth_selected |> select(fam, id, Sex, Dx)
    covar_merged |>
        write.table(file=out_cov, sep="\t", row.names=TRUE,
                    col.names=FALSE, quote=FALSE)
    qcovar_merged <- meth_selected |> select(fam, id, Age)
    qcovar_merged |>
        write.table(file=out_qcov, sep="\t", row.names=TRUE,
                    col.names=FALSE, quote=FALSE)
  return(list(covar_merged=covar_merged, qcovar_merged=qcovar_merged))
}

## Main
                                        # load data
load(here("inputs/wgbs-data/caudate", paste0("Caudate_chr",chr,"_BSobj.rda")))
output_path <- here("testing/caudate/_m", paste0("chr_",chr))

                                        # create output directory if it  
                                        # doesn't exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

                                        # change file path for raw data
raw_assays  <- here("inputs/wgbs-data/caudate/raw/CpGassays.h5")
BSobj       <- change_file_path(BSobj,raw_assays)

                                        # keep only adult AA
pheno_file_path <- here("inputs/phenotypes/merged/_m/merged_phenotypes.csv")
BSobj           <- filter_pheno(BSobj,pheno_file_path)

                                        # exclude low coverage sites
BSobj <- exclude_low_cov(BSobj)

                                        # calculate sd & mean of DNAm
##dir.create()
stats <- DNAm_stats(BSobj,file.path(output_path, "stats.rda"))

                                        # extract VMRs
v   <- data.frame(chr=chr,start=start(BSobj),sd=stats$sds)
vmr <- get_vmr(v)

                                        # read in FID, IID from sample file
psam_file <- here("inputs/genotypes/TOPMed_LIBD.AA.psam")
samples   <- extract_fid_iid(psam_file)

                                        # merge methylation values with FID and
                                        # IID and write to .phen file
meth_merged <- write_meth_to_phen(BSobj, stats$M, samples, output)

# write covariate files
covars <- write_covar(pheno_file_path,meth_merged,output_path)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
