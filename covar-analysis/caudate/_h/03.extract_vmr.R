#### Extract VMRs ####

suppressPackageStartupMessages({
    library('data.table')
    library(here)
})

args <- commandArgs(trailingOnly = TRUE)
chr  <- args[1]

# Function
get_vmr <- function(v, out_vmr) {
  sdCut  <- quantile(v[, 3], prob = 0.99, na.rm = TRUE)
  vmrs   <- c()
  v2     <- v[v$chr==chr, ]
  v2     <- v2[order(v2$start), ]
  isHigh <- rep(0, nrow(v2))
  isHigh[v2$sd > sdCut] <- 1
  vmrs0  <- bsseq:::regionFinder3(isHigh, as.character(v2$chr), 
                                  v2$start, maxGap = 1000)$up
  vmr    <- vmrs0[vmrs0$n > 5,1:3]
  write.table(vmr, file=out_vmr, col.names=F,
              row.names=F, sep="\t", quote=F)
  return(vmr)
}

# Main

                                        # create output dir if doesn't exist
output_path <- here("covar-analysis/caudate/_m/vmr", paste0("chr_", chr))

if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

                                        # load sd of pc residuals
var_file <- here("covar-analysis", "caudate", "_m", "pca", 
                 paste0("chr_", chr), "res_var_all.tsv")
res_var <- fread(var_file, select = 1:3, header = TRUE)
colnames(res_var) <- c("chr", "start", "sd")

                                        # extract VMRs
vmr <- get_vmr(res_var, file.path(output_path, "vmr.bed"))

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()