#### Calculate mean and SD of methylation values ####

library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')
library(scales)
here::i_am("testing/_h/01.sd_test.R")
library(here)
library(dplyr)
library(genio)
library(tibble)
library(plinkr)

### Adapting sd.R 01_pca.R and 03_vmr.R for one chromosome
setwd("/projects/p32505/projects/dna-methylation-heritability/testing/_h")
load(here("inputs/wgbs-data/caudate/Caudate_chr1_BSobj.rda"))

# Change file path for raw data
BSobj@assays@data@listData$M@seed@seed@filepath <- here("inputs/wgbs-data/caudate/raw/CpGassays.h5")
BSobj@assays@data@listData$Cov@seed@seed@filepath <- here("inputs/wgbs-data/caudate/raw/CpGassays.h5")
BSobj@assays@data@listData$coef@seed@seed@filepath <- here("inputs/wgbs-data/caudate/raw/CpGassays.h5")

#BSobj = BSobj[seqnames(BSobj) == "chr1",] #keep chr1

# keep only control, adult AA
pheno <- here("inputs/phenotypes/merged/_m/merged_phenotypes.csv")
ances <- read.csv(pheno,header=T)
id <- intersect(ances$BrNum[ances$Race == "AA" & ances$Dx == "Control" & ances$Age >= 17],colData(BSobj)$brnum)
BSobj <- BSobj[,is.element(colData(BSobj)$brnum,id)]

# exlcude low coverage sites
cov=getCoverage(BSobj)
n <- length(colData(BSobj)$brnum)
keep <- which(rowSums2(cov >= 5) >= n * 0.8) 
BSobj <- BSobj[keep,]

# calculate sd & mean of DNAm
M=as.matrix(getMeth(BSobj))
sds <- rowSds(M)
means <- rowMeans2(M)
save(sds,means,BSobj,file=paste0("./","chr1_stats",".rda"))

# read in FID, IID from sample file 
samples <- read_plink2_psam_file("/projects/p32505/projects/dna-methylation-heritability/inputs/genotypes/TOPMed_LIBD.AA.psam")
samples <- samples[,1:2]

# add sample IDs to methylation matrix
sample_ids <- colData(BSobj)$brnum #use this or id?
M_transpose <- t(M)
meth_df <- data.frame(FID = sample_ids, M_transpose)

# merge methylation values with FID and IID
meth_merged <- meth_df %>% 
  inner_join(samples, by = "FID") %>% 
  arrange(match(FID, samples$FID))
meth_merged <- meth_merged %>% 
  select(FID, IID, everything())

#Write methylation values to .phen file
colnames(meth_merged)[1:3] <- c("fam", "id", "pheno")
write_phen("methylation.phen", meth_merged)