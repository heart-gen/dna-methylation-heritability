##### Summarizes heritability statistics at different filtering thresholds #####
suppressPackageStartupMessages({
    library(here)
    library(qvalue)
    library(dplyr)
    library(tidyverse)
})

# Define brain regions and p-value thresholds
tissues         <- c("hippocampus", "caudate", "dlpfc") 
pval_thresholds <- c(0.05, 0.1, 0.25)

# Loop over tissues
results <- map(tissues, function(tissue) {
  
    # Read in vmrs
    vmr <- read.table(here("heritability/gcta/",
                            paste0(tissue, "/_m/summary/greml_summary.tsv")),
                        sep = "\t", header = TRUE)
    n_total_vmrs <- nrow(vmr)
    
    # FDR correction
    vmr$p_adjusted_fdr <- p.adjust(vmr$Pval_Variance, method = "fdr")
    vmr$q_val_fdr      <- qvalue(vmr$Pval_Variance, lambda = seq(0, 0.5, 0.01))$qvalues
    
    # Get top 10% heritable sites
    top_10        <- quantile(vmr$Sum.of.V.G._Vp_Variance, 0.90, na.rm = TRUE)
    vmr_heritable <- filter(vmr, Sum.of.V.G._Vp_Variance >= top_10)
    
    # Summary stats for heritable sites at each p-value threshold
    stat_tbl <- map_dfr(pval_thresholds, function(thresh) {
        df_filt <- filter(vmr_heritable, p_adjusted_fdr < thresh)
        tibble(
        tissue = tissue,
        n_total_vmrs = n_total_vmrs,
        top_10_thresh = as.character(top_10),
        p_thresh = as.character(thresh),
        n_significant = nrow(df_filt),
        mean_h2 = mean(df_filt$Sum.of.V.G._Vp_Variance, na.rm = TRUE),
        median_h2 = median(df_filt$Sum.of.V.G._Vp_Variance, na.rm = TRUE)
        )
    })
    
    # Summary stats for all VMRs
    all_stats <- tibble(
        tissue = tissue,
        n_total_vmrs = n_total_vmrs,
        top_10_thresh = "all",
        p_thresh = "all",
        n_significant = NA,
        mean_h2 = mean(vmr$Sum.of.V.G._Vp_Variance, na.rm = TRUE),
        median_h2 = median(vmr$Sum.of.V.G._Vp_Variance, na.rm = TRUE)
    )
    bind_rows(all_stats, stat_tbl)
})

# Combine results and write to file
summary_df <- bind_rows(results)
write_tsv(summary_df, file = here("heritability/gcta/tissue_comparison/summary/_m", "summary_h2_filtering.tsv"))
