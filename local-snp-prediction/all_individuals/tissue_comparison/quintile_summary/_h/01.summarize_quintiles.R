## This script summarizes h2 values and concordance for quintiles of matched VMRs in AA and EA

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

N_QUINTILES <- 5

## Function
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

  quintile_h2_summary_AA <- df %>%
    group_by(h2_quintile_AA) %>%
    summarise(
        n_vmrs = n(),
        h2_min = min(h2_unscaled_AA, na.rm = TRUE),
        h2_max = max(h2_unscaled_AA, na.rm = TRUE),
        h2_mean = mean(h2_unscaled_AA, na.rm = TRUE),
        h2_sd = sd(h2_unscaled_AA, na.rm = TRUE)
    ) %>% rename(quintile = h2_quintile_AA)

  write.csv(quintile_h2_summary_AA, 
            file = file.path(out_path, 
                             paste0("quintile_h2_summary_AA_", tissue, ".csv")), 
            row.names = FALSE)

  quintile_h2_summary_EA <- df %>%
    group_by(h2_quintile_EA) %>%
    summarise(
        n_vmrs = n(),
        h2_min = min(h2_unscaled_EA, na.rm = TRUE),
        h2_max = max(h2_unscaled_EA, na.rm = TRUE),
        h2_mean = mean(h2_unscaled_AA, na.rm = TRUE),
        h2_sd = sd(h2_unscaled_AA, na.rm = TRUE)
    ) %>% rename(quintile = h2_quintile_EA)
  
  write.csv(quintile_h2_summary_EA, 
            file = file.path(out_path, 
                             paste0("quintile_h2_summary_EA_", tissue, ".csv")), 
            row.names = FALSE)

  quintile_match_summary <- df %>%
    count(h2_quintile_AA, h2_quintile_EA) %>%
    tidyr::complete(h2_quintile_AA, h2_quintile_EA, fill = list(n = 0))
  
  matched_count <- sum(df$h2_quintile_AA == df$h2_quintile_EA, na.rm = TRUE)
  percent_matched <- (matched_count / nrow(df)) * 100
  
  print(paste0(matched_count, " VMRs matched across quintiles (", 
               round(percent_matched, 2), "%)"))
  
  write.csv(quintile_match_summary, 
            file = file.path(out_path, 
                             paste0("quintile_match_summary_AA_EA_", tissue, ".csv")), 
            row.names = FALSE)
}

## Main

tissues <- c("caudate", "hippocampus", "dlpfc")

out_path <- here("heritability/elastic_net_model/all_individuals/tissue_comparison/quintile_summary/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

for (tissue in tissues) {
  # Read in summary table
  enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net_matched_r2_0.3.tsv"))
  enet <- read.table(enet_file, sep = "\t", header = TRUE)
  
  # Quantify matched VMRs across quintiles
  quintile_summary(enet, tissue, out_path)
}

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()