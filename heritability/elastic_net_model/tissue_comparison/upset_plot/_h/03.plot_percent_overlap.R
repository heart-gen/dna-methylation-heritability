## This script plots the distribution of percent overlap values

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(tibble)
  library(tidyr)
  library(ggpubr)
})

## Function
get_overlaps <- function(flag, tissue1, tissue2, h2_cat){
  
  overlap_fn <- paste0(flag, "/percent_overlap/", tissue1, "_", tissue2, "_overlap_", h2_cat, ".tsv")
  pct_overlap <- fread(overlap_fn, sep = "\t", header = TRUE)
  
  return(pct_overlap)
}

save_plot <- function(p, fn, w, h, dpi){
  for(ext in c(".pdf", ".png")){
    ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h, dpi=dpi)
  }
}

plot_dist <- function(pct_overlap, tissue1, tissue2, h2_cat, flag){
  pct_long <- pct_overlap %>%
    pivot_longer(cols = c("pctA", "pctB", "reciprocal"),
                 names_to = "option", values_to = "percent")
  
  threshold_pct <- (as.numeric(sub(".*_", "", flag)) * 100)
  
  p <- gghistogram(pct_long, x = "percent",
                   add_density = TRUE, rug = TRUE,
                   add = "median", fill = "option",
                   color = "option") +
    facet_wrap(~option) + 
    theme_pubr(base_size = 15, border = TRUE) +
    labs(
      title = paste("Percent overlap distribution for", tissue1, "vs.", tissue2),
      subtitle = paste(h2_cat, "VMRs"),
      x = "Percent overlap (%)"
    ) +
    geom_vline(xintercept = threshold_pct, linetype = "dashed", color = "black") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
}

## Main
out_path <- here("heritability/elastic_net_model/tissue_comparison/upset_plot/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

option_flags <- c("F_0.25", "F_0.5", "F_0.75", "f_0.25", "f_0.5", "f_0.75")
h2_categories <- c("heritable", "non-heritable", "low_prediction", "all")

for (flag in option_flags) {
  for (h2_cat in c("heritable", "non-heritable", "low_prediction")){
    for (tissue1 in c("caudate", "hippocampus")){
      for (tissue2 in c("hippocampus", "dlpfc")){
        if (tissue1 != tissue2){
          pct_overlap <- get_overlaps(flag, tissue1, tissue2, h2_cat)
          out_dist <- file.path(out_path, flag, "percent_overlap",
                                paste0("VMR_distribution_", h2_cat))
          p <- plot_dist(pct_overlap, tissue1, tissue2, h2_cat, flag)
          save_plot(p, out_dist, 8, 6, 300)
        }
      }
    }
  }
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()