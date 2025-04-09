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
chr <- args[1]

here::i_am("testing/caudate/_h/01.sd_test.R")
setwd(here("testing/caudate/_h"))
load(here(paste0("inputs/wgbs-data/caudate/Caudate_chr",chr,"_BSobj.rda")))
output_path <- here(paste0("testing/caudate/_m/chr_",chr))

# change file path for raw data
change_file_path <- function(BSobj,raw_assays) {
  var <- c("M", "Cov", "coef")
  for (assay in var) {
    BSobj@assays@data@listData[[assay]]@seed@seed@filepath <- raw_assays
  }
  return(BSobj)
}
raw_assays <- here("inputs/wgbs-data/caudate/raw/CpGassays.h5")
BSobj <- change_file_path(BSobj,raw_assays)

# keep only adult AA
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
pheno_file_path <- here("inputs/phenotypes/merged/_m/merged_phenotypes.csv")
BSobj <- filter_pheno(BSobj,pheno_file_path)

# exclude low coverage sites
exclude_low_cov <- function(BSobj) {
  cov=getCoverage(BSobj)
  n <- length(colData(BSobj)$brnum)
  keep <- which(rowSums2(cov >= 5) >= n * 0.8) 
  BSobj <- BSobj[keep,]
  return(BSobj)
}
BSobj <- exclude_low_cov(BSobj)

# calculate sd & mean of DNAm
DNAm_stats <- function(BSobj,output_stats) {
  M=as.matrix(getMeth(BSobj))
  sds <- rowSds(M)
  means <- rowMeans2(M)
  save(sds,means,BSobj,file=output_stats)
  return(list(M = M, sds = sds, means = means))
}
stats <- DNAm_stats(BSobj,file.path(output_path, "stats.rda"))

# read in FID, IID from sample file 
extract_fid_iid <- function(psam_file) {
  samples <- read_plink2_psam_file(psam_file)
  samples <- samples[,1:2]
  return(samples)
}
psam_file <- here("inputs/genotypes/TOPMed_LIBD.AA.psam")
samples <- extract_fid_iid(psam_file)

# merge methylation values with FID and IID and write to .phen file
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
meth_merged <- write_meth_to_phen(BSobj,stats$M,samples,output)

# write covariate files
write_covar <- function(pheno_file_path,meth_merged,output_path) {
  pheno <- read.csv(pheno_file_path,header=T)
  id <- intersect(
    pheno$BrNum[
      pheno$Race == "AA" & 
      pheno$Age >= 17 & 
      pheno$Region == "Caudate"
    ], 
    colData(BSobj)$brnum
  )
  covar <- pheno %>%
    filter(Region == "Caudate") %>%
    select(BrNum, Sex, Dx)
  qcovar <- pheno %>%
    filter(Region == "Caudate") %>%
    select(BrNum, Age)
  covar_filtered <- covar %>%
    filter(BrNum %in% id)
  qcovar_filtered <- qcovar %>%
    filter(BrNum %in% id)
  meth_selected <- meth_merged %>%
    select(fam, id)
  covar_merged <- meth_selected %>%
    inner_join(covar, by = c("fam" = "BrNum")) %>%
    arrange(match(fam, meth_selected$fam))
  qcovar_merged <- meth_selected %>%
    inner_join(qcovar, by = c("fam" = "BrNum")) %>%
    arrange(match(fam, meth_selected$fam))
  write.table(
    covar_merged, 
    file = file.path(output_path, "TOPMed_LIBD.AA.covar"), 
    sep = "\t", 
    row.names = FALSE, 
    col.names = FALSE, 
    quote = FALSE
  )
  write.table(
    qcovar_merged, 
    file = file.path(output_path, "TOPMed_LIBD.AA.qcovar"), 
    sep = "\t", 
    row.names = FALSE, 
    col.names = FALSE, 
    quote = FALSE
  )
  return(list(covar_merged=covar_merged, qcovar_merged=qcovar_merged))
}
covars <- write_covar(pheno_file_path,meth_merged,output_path)

#### Reproducibility information
if(sge_id == 1){
  Sys.time()
  proc.time()
  options(width = 120)
  sessioninfo::session_info()
}