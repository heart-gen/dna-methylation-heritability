#### Summarize top 1% variable cpg sites per brain region ####

suppressPackageStartupMessages({
  library('data.table')
  library(here)
})

## Main

chromosomes <- c(1:22, "X", "Y")

total_top1 <- 0
results <- list()

for (chr in chromosomes) {
  
  cat("Processing chr", chr, "\n")
  
  var_file <- here("vmr-analysis", "caudate", "_m", "pca",
                   paste0("chr_", chr), "res_var_all.tsv")
  
  if (!file.exists(var_file)) {
    cat("Residual variance file not found, skipping\n")
    next
  }
  
  v <- fread(var_file, select = 1:3, header = TRUE)
  colnames(v) <- c("chr", "start", "sd")
  
  sdCut <- quantile(v$sd, probs = 0.99, na.rm = TRUE)
  n_top1 <- sum(v$sd > sdCut, na.rm = TRUE)
  
  total_top1 <- total_top1 + n_top1
  results[[as.character(chr)]] <- data.table(
    tissue = "Caudate", chr = chr, total_cpg = nrow(v), 
    top1_cpg = n_top1, sd_cutoff = sdCut
  )
}

results_dt <- rbindlist(results)
print(results_dt)

out_cpg <- here("vmr-analysis/caudate/_m/cpg")
write.table(results_dt, file=file.path(out_cpg, "top1_cpg.tsv"), 
            col.names=T, row.names=F, sep="\t", quote=F)

cat("Total top 1% CpGs across all chr:", total_top1)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()