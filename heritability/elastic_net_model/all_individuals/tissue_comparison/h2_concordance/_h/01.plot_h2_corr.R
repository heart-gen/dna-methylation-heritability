#### Plot h2 concordance of matched VMRs in AA and EA ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(tidyverse)
  library(ggpubr)
})

N_QUINTILES <- 5

## Function
filter_sites <- function(enet) {
  vmr <- na.omit(enet)
  vmr <- vmr %>%
    mutate(h2_category = case_when(
      r_squared_cv <= 0.3 ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.3 ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.3 ~ "Heritable"
    ),
    h2_category = factor(h2_category, levels = c("Heritable", 
                                                 "Non-heritable", 
                                                 "Low prediction"))
    )
  return(vmr)
}

quintile_summary <- function(df, tissue, out_path){
  df <- df %>%
    filter(h2_category %in% c("Heritable", "Non-heritable"))
  
  breaks_AA <- quantile(df$h2_unscaled_AA,
                        probs = seq(0, 1, 1 / N_QUINTILES),
                        na.rm = TRUE) |>
    unique()
  n_bins_AA <- length(breaks_AA) - 1
  
  breaks_EA <- quantile(df$h2_unscaled_EA,
                        probs = seq(0, 1, 1 / N_QUINTILES),
                        na.rm = TRUE) |>
    unique()
  n_bins_EA <- length(breaks_EA) - 1
  
  df <- df |>
    mutate(h2_quintile_AA = cut(h2_unscaled_AA,
                                breaks = breaks_AA,
                                labels = paste0("Q", seq_len(n_bins_AA)),
                                include.lowest = TRUE)) |>
    mutate(h2_quintile_EA = cut(h2_unscaled_EA,
                                breaks = breaks_EA,
                                labels = paste0("Q", seq_len(n_bins_EA)),
                                include.lowest = TRUE)) |>
    filter(!is.na(h2_quintile_AA), !is.na(h2_quintile_EA))
  
  quintile_summary <- df %>%
    count(h2_quintile_AA, h2_quintile_EA) %>%
    tidyr::complete(h2_quintile_AA, h2_quintile_EA, fill = list(n = 0))
  
  matched_count <- sum(df$h2_quintile_AA == df$h2_quintile_EA, na.rm = TRUE)
  percent_matched <- (matched_count / nrow(df)) * 100
  
  print(paste0(matched_count, " VMRs matched across quintiles (", 
               round(percent_matched, 2), "%)"))
  
  write.csv(quintile_summary, 
            file = file.path(out_path, 
                             paste0("quintile_match_summary_AA_EA_", tissue, ".csv")), 
            row.names = FALSE)
}

spearman_corr <- function(h2_df, h2_cat, tissue, out_path){
  spearman <- h2_df %>% 
    filter(h2_category == h2_cat) %>%
    summarise(
      spearman_rho = cor.test(h2_unscaled_AA, h2_unscaled_EA, method = "spearman")$estimate,
      spearman_p_value = cor.test(h2_unscaled_AA, h2_unscaled_EA, method = "spearman")$p.value,
      n = n()
    )
  print(spearman)
  write.csv(spearman, 
            file = file.path(out_path, 
                             paste0(gsub(" ", "_", tolower(h2_cat)), "_h2_corr_AA_EA_", 
                                    tissue, ".csv")), 
            row.names = FALSE)
}

plot_corr <- function(h2_df, h2_cat, tissue, output_path){
  h2_df <- h2_df %>%
    filter(h2_category == h2_cat)
  
  xlab = "h2_BA"
  ylab = "h2_WA"
  
  fn = file.path(output_path, paste(gsub(" ", "_", tolower(h2_cat)), "h2_corr_AA_EA",
                                    tissue, sep="_"))
  
  pp <- ggscatter(h2_df, x = "h2_unscaled_AA", y = "h2_unscaled_EA", add = "reg.line", 
                  size = 1, xlab = xlab, ylab = ylab, panel.labs.font=list(face="bold"),
                  add.params=list(color = "blue", fill = "lightgray"),
                  conf.int = TRUE, cor.coef = TRUE, cor.coef.size = 3,
                  cor.method = "spearman", cor.coeff.args=list(label.sep="\n"), ncol = 4) +
    theme_pubr(base_size=18) + 
    labs(title = paste("BA vs WA h2 correlation:", h2_cat))
  save_plot(pp, fn, 6, 6)
}

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}

tissues <- c("caudate", "hippocampus", "dlpfc")

out_path <- here("heritability/elastic_net_model/all_individuals/tissue_comparison/h2_concordance/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

for (tissue in tissues) {
  # Read in summary table
  enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net_AA.tsv"))
  enet_AA <- read.table(enet_file, sep = "\t", header = TRUE)
  
  enet_file_EA <- here("heritability/elastic_net_model/all_individuals/", 
                       paste0(tissue, "/_m/", tissue, "_summary_elastic-net_EA.tsv"))
  enet_EA <- read.table(enet_file_EA, sep = "\t", header = TRUE)
  
  vmr_AA <- filter_sites(enet_AA)
  vmr_EA <- filter_sites(enet_EA)
  
  # Get matched VMRs by h2 cat
  vmr_combined <- inner_join(vmr_AA, vmr_EA, by = c("chrom", "start", "end"),
                             suffix = c("_AA", "_EA")) %>%
    filter(h2_category_AA == h2_category_EA) %>%
    rename(h2_category = h2_category_AA) %>%
    select(-c(h2_category_EA, race_AA, race_EA, region_EA))
  
  vmr_combined <- vmr_combined %>%
    mutate(feature_id = paste(chrom, start, end, sep = "_"))
  
  # Quantify matched VMRs across quintiles
  quintile_summary(vmr_combined, tissue, out_path)
  
  h2_cats <- c("Heritable", "Non-heritable", "Low prediction")
  
  for (h2_cat in h2_cats) {
    
    # Spearman correlation
    spearman_corr(vmr_combined, h2_cat, tissue, out_path)
    
    # Plot scatter for matched VMRs
    plot_corr(vmr_combined, h2_cat, tissue, out_path)
    
  }
}

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()