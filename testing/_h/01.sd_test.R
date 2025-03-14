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

#Write methylation values to .phen file
meth=as_tibble(getMeth(BSobj))
average_meth <- tibble(pheno = rowMeans(meth)) 
average_meth <- average_meth %>%
  mutate(
    fam = rep(NA, nrow(average_meth)),
    id = rep(NA, nrow(average_meth)) 
  ) %>%
  select(fam, id, everything())
write_phen("methylation.phen", meth)