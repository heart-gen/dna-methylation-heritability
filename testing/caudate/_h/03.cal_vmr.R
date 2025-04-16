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

args      <- commandArgs(trailingOnly = TRUE)
chr_num   <- args[1]
chr       <- paste0("chr", chr_num)
start_pos <- as.integer(args[2])
end_pos   <- as.integer(args[3])

#here::i_am("testing/caudate/_h/03.cal_vmr.R")
## Function
calc_vmr_meth <- function(BSobj, chr, start_pos, end_pos) {
  reg      <- GRanges(seqnames = chr, 
                      ranges = IRanges(start = start_pos, end = end_pos))
  meth_reg <- getMeth(BSobj, regions = reg, what="perRegion")
  return(meth_reg)
}

## Main
# load raw DNAm
load(here("testing/caudate/_m", paste0("chr_", chr_num, "/stats.rda")))
output_path <- here("testing/caudate/_m", paste0("chr_", chr_num))

# extract methylation values for each VMR
meth_reg <- calc_vmr_meth(BSobj, chr, start_pos, end_pos)
  
# read in FID, IID from sample file
psam_file <- here("inputs/genotypes/TOPMed_LIBD.AA.psam")
samples   <- extract_fid_iid(psam_file)

# merge methylation values with FID and
# IID and write to .phen file
out_phen    <- file.path(output_path, 
                         paste0(start_pos, "_", end_pos, "_meth.phen"))
meth_merged <- write_meth_to_phen(BSobj, meth_reg, samples, out_phen)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
