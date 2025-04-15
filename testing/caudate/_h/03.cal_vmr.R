#### Calculate methylation values for variably methylated regions ####

library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')
library(scales)
library(GenomicRanges)
library(dplyr)
library(genio)
library(plinkr)
library(here)
source(here("testing/caudate/_h/01.sd_test.R"))

args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
start_pos <- args[2]
end_pos <- args[3]

#here::i_am("testing/caudate/_h/03.cal_vmr.R")
output_path <- here(paste0("testing/caudate/_m/chr_",chr))

# load raw DNAm
load(here(paste0("testing/caudate/_m/chr_",chr,"/stats.rda")))

# extract methylation values for each VMR
calc_vmr_meth <- function(BSobj,chr,start_pos,end_pos) {
  reg      <- GRanges(seqnames = chr, 
                      ranges = IRanges(start = start_pos, end = end_pos))
  meth_reg <- getMeth(BSobj, regions = reg, what="perRegion")
  return(meth_reg)
}

# main
meth_reg <- calc_vmr_meth(BSobj,chr,start_pos,end_pos)
  
# read in FID, IID from sample file
psam_file <- here("inputs/genotypes/TOPMed_LIBD.AA.psam")
samples   <- extract_fid_iid(psam_file)

# merge methylation values with FID and
# IID and write to .phen file
out_phen    <- file.path(output_path, 
                         paste0(start_pos,"_",end_pos,"_meth.phen"))
meth_merged <- write_meth_to_phen(BSobj, meth_reg, samples, out_phen)


# read in FID, IID from sample file 
#samples <- read_plink2_psam_file(here("inputs/inputs/genotypes/TOPMed_LIBD.AA.psam"))
#samples <- samples[,1:2]

# add sample IDs to methylation matrix
#sample_ids <- colData(BSobj)$brnum #use this or id?
#meth_reg_transpose <- t(meth_reg)
#meth_reg_df <- data.frame(FID = sample_ids, meth_reg_transpose)

# merge methylation values with FID and IID
#meth_reg_merged <- meth_reg_df %>% 
#  inner_join(samples, by = "FID") %>% 
#  arrange(match(FID, samples$FID))
#meth_reg_merged <- meth_reg_merged %>% 
#  select(FID, IID, everything())

# write methylation values to .phen file
#colnames(meth_reg_merged)[1:3] <- c("fam", "id", "pheno")
#write_phen(file=file.path(output_path,"VMR1_meth.phen"), meth_reg_merged)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
