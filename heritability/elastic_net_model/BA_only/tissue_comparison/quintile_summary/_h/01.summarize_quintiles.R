## This script summarizes h2 values for quintiles in BA only cohort

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
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
  
  breaks <- quantile(df$h2_unscaled,
                        probs = seq(0, 1, 1 / N_QUINTILES),
                        na.rm = TRUE) |>
    unique()
  n_bins <- length(breaks) - 1
  
  df <- df |>
    mutate(h2_quintile = cut(h2_unscaled,
                                breaks = breaks,
                                labels = paste0("Q", seq_len(n_bins)),
                                include.lowest = TRUE)) |>
    filter(!is.na(h2_quintile))

  quintile_h2_summary <- df %>%
    group_by(h2_quintile) %>%
    summarise(
        n_vmrs = n(),
        h2_min = min(h2_unscaled, na.rm = TRUE),
        h2_max = max(h2_unscaled, na.rm = TRUE),
        h2_mean = mean(h2_unscaled, na.rm = TRUE),
        h2_sd = sd(h2_unscaled, na.rm = TRUE)
    ) %>% rename(quintile = h2_quintile)

  write.csv(quintile_h2_summary, 
            file = file.path(out_path, 
                             paste0("quintile_h2_summary_", tissue, ".csv")), 
            row.names = FALSE)
}

## Main

tissues <- c("caudate", "hippocampus", "dlpfc")

out_path <- here("heritability/elastic_net_model/BA_only/tissue_comparison/quintile_summary/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

for (tissue in tissues) {
  # Read in summary table
  enet_file <- here("heritability/elastic_net_model/BA_only/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
  enet <- read.table(enet_file, sep = "\t", header = TRUE)

  vmr <- filter_sites(enet)
  
  # Quantify matched VMRs across quintiles
  quintile_summary(vmr, tissue, out_path)
}

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()