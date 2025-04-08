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

### Adapting sd.R 01_pca.R and 03_vmr.R for one chromosome
here::i_am("testing/caudate/_h/01.sd_test.R")
setwd(here("testing/caudate/_h"))
load(here("inputs/wgbs-data/caudate/Caudate_chr1_BSobj.rda"))
output <- here("testing/caudate/_m/chr_1")

# Change file path for raw data
BSobj@assays@data@listData$M@seed@seed@filepath <- here("inputs/wgbs-data/caudate/raw/CpGassays.h5")
BSobj@assays@data@listData$Cov@seed@seed@filepath <- here("inputs/wgbs-data/caudate/raw/CpGassays.h5")
BSobj@assays@data@listData$coef@seed@seed@filepath <- here("inputs/wgbs-data/caudate/raw/CpGassays.h5")

# keep only control, adult AA
pheno <- here("inputs/phenotypes/merged/_m/merged_phenotypes.csv")
ances <- read.csv(pheno,header=T)
id <- intersect(ances$BrNum[ances$Race == "AA" & ances$Dx == "Control" & ances$Age >= 17 & ances$Region == "Caudate"],colData(BSobj)$brnum)
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
save(sds,means,BSobj,file=file.path(output, "chr_1_stats.rda"))

# read in FID, IID from sample file 
psam <- here("inputs/genotypes/TOPMed_LIBD.AA.psam")
samples <- read_plink2_psam_file(psam)
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

# write methylation values to .phen file
colnames(meth_merged)[1:3] <- c("fam", "id", "pheno")
write_phen(file=file.path(output,"chr_1_cpg_meth.phen"), meth_merged)

# add sample IDs to covariate file
covar <- ances %>%
  filter(Region == "Caudate") %>%
  select(BrNum, Sex)
qcovar <- ances %>%
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
write.table(covar_merged, file=file.path(output,"chr_1.covar"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(qcovar_merged, file=file.path(output,"chr_1.qcovar"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
