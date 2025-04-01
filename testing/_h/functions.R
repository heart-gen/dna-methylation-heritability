#### VMR Test ####

library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')
library(scales)
library(sessioninfo)
library(here)

### Adapting sd.R 01_pca.R and 03_vmr.R for one chromosome
load_bsobj <- function(file_name) {
  load(file_name)     # Load the .rda file
}

load("Caudate_chr1_BSobj.rda")
#load_bsobj("Caudate_chr1_BSobj.rda")
print(exists("BSobj"))
