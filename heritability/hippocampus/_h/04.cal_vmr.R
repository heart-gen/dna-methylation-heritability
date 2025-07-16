#### Calculate methylation values for variably methylated regions ####

suppressPackageStartupMessages({
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
})

                                        # get loci from command line
args      <- commandArgs(trailingOnly = TRUE)
chr_num   <- args[1]
chr       <- paste0("chr", chr_num)
start_pos <- as.integer(args[2])
end_pos   <- as.integer(args[3])

## Function
calc_vmr_meth <- function(BSobj, chr, start_pos, end_pos) {
  reg      <- GRanges(seqnames = chr, 
                      ranges = IRanges(start = start_pos, end = end_pos))
  meth_reg <- getMeth(BSobj, regions = reg, what="perRegion")
  return(meth_reg)
}

extract_fid_iid <- function(psam_file) {
    samples <- read_plink2_psam_file(psam_file)
    samples <- samples[, 1:2]
    return(samples)
}

write_meth_to_phen <- function(BSobj, M, samples, out_phen) {

                                        # add sample IDs to methylation matrix
    sample_ids <- colData(BSobj)$brnum
    meth_df    <- data.frame(FID = sample_ids, t(M))

                                        # merge methylation data and 
                                        # sample ids by FID
    meth_merged <- meth_df %>%
        inner_join(samples, by = "FID") %>%
        arrange(match(FID, samples$FID)) %>%
        select(FID, IID, everything())
    colnames(meth_merged)[1:3] <- c("fam", "id", "pheno")

                                        # write methylation values to .phen file
    write_phen(file=out_phen, meth_merged)
    return(meth_merged)
}

## Main

                                        # load raw DNAm
load(here("heritability/hippocampus/_m/cpg", paste0("chr_", chr_num, "/stats.rda")))
output_path <- here("heritability/hippocampus/_m/vmr", paste0("chr_", chr_num))

                                        # calculate DNAm for each VMR
meth_reg <- calc_vmr_meth(BSobj, chr, start_pos, end_pos)

                                        # read FID, IID from sample file
psam_file <- here("inputs/genotypes/TOPMed_LIBD.AA.psam")
samples   <- extract_fid_iid(psam_file)

                                        # merge DNAm with FID, IID
out_phen    <- file.path(output_path, 
                         paste0(start_pos, "_", end_pos, "_meth.phen"))

                                        # write to .phen file
meth_merged <- write_meth_to_phen(BSobj, meth_reg, samples, out_phen)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
