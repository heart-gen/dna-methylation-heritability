#### Extract variably methylated regions ####

library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')
library(scales)
library(GenomicRanges)
library(dplyr)
library(genio)
library(plinkr)

args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
start_pos <- args[2]
end_pos <- args[3]

here::i_am("testing/caudate/_h/01.sd_test.R")
output_path <- here(paste0("testing/caudate/_m/chr_",chr))

# load sd of raw DNAm
load(here(paste0("testing/caudate/_m/chr_",chr,"/stats.rda")))
v <- data.frame(chr=chr,start=start(BSobj),sd=sds)

# sd cutoff
sdCut <- quantile(v[,3], prob = 0.99, na.rm = TRUE)

# get vmr
vmrs <- c()
v2 <- v[v$chr==chr,]
v2 <- v2[order(v2$start),]
isHigh <- rep(0, nrow(v2))
isHigh[v2$sd > sdCut] <- 1
vmrs0 <- bsseq:::regionFinder3(isHigh, as.character(v2$chr), v2$start, maxGap = 1000)$up
vmr <- vmrs0[vmrs0$n > 5,1:3]
write.table(vmr,file=file.path(output_path,"vmr.bed"),col.names=F,row.names=F,sep="\t",quote=F)

# extract methylation values for each VMR
reg <- GRanges(seqnames = chr, ranges = IRanges(start = start_pos, end = end_pos))
meth_reg <- getMeth(BSobj, regions = reg, what="perRegion")

# read in FID, IID from sample file 
samples <- read_plink2_psam_file(here("inputs/inputs/genotypes/TOPMed_LIBD.AA.psam"))
samples <- samples[,1:2]

# add sample IDs to methylation matrix
sample_ids <- colData(BSobj)$brnum #use this or id?
meth_reg_transpose <- t(meth_reg)
meth_reg_df <- data.frame(FID = sample_ids, meth_reg_transpose)

# merge methylation values with FID and IID
meth_reg_merged <- meth_reg_df %>% 
  inner_join(samples, by = "FID") %>% 
  arrange(match(FID, samples$FID))
meth_reg_merged <- meth_reg_merged %>% 
  select(FID, IID, everything())

# write methylation values to .phen file
colnames(meth_reg_merged)[1:3] <- c("fam", "id", "pheno")
write_phen(file=file.path(output_path,"VMR1_meth.phen"), meth_reg_merged)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
