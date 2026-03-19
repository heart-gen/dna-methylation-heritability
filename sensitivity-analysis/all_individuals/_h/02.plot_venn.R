# Summarize heritability classifications using stacked barplot

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggvenn)
  library(tidyverse)
  library(ggpubr)
})

## Function
save_plot <- function(p, fn, w, h, dpi){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn,ext), plot=p, width=w, height=h, dpi=dpi)
  }
}

## Main 
# Create output directory
out_path <- here("sensitivity-analysis/all_individuals/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

tissues <- c("caudate", "hippocampus", "dlpfc")
populations <- c("AA", "EA")

for (pop in populations) {
  
  plot_list <- list()
  
  for (tissue in tissues) {
    # Read in summary table
    enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                      paste0(tissue, "/_m/", tissue, "_summary_elastic-net_", pop, ".tsv"))
    enet <- read.table(enet_file, sep = "\t", header = TRUE)
    
    # Get summary for each r2 threshold
    enet_high <- enet %>%
      mutate(h2_category = case_when(
        r_squared_cv <= 0.75 ~ "Low prediction",
        h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
        h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
      ),
      r2_threshold = 0.75)
    
    enet_high <- enet_high %>%
      mutate(h2_category = factor(h2_category,
                                  levels = c("Heritable", "Non-heritable", "Low prediction")),
             region = recode(region, "caudate" = "Caudate"))
    
    enet_low <- enet %>%
      mutate(h2_category = case_when(
        r_squared_cv <= 0.3 ~ "Low prediction",
        h2_unscaled < 0.1 & r_squared_cv > 0.3 ~ "Non-heritable",
        h2_unscaled >= 0.1 & r_squared_cv > 0.3 ~ "Heritable"
      ),
      r2_threshold = 0.3)
    
    enet_low <- enet_low %>%
      mutate(h2_category = factor(h2_category,
                                  levels = c("Heritable", "Non-heritable", "Low prediction")),
             region = recode(region, "caudate" = "Caudate"))
    
    all_summary <- bind_rows(enet_high, enet_low)
    
    all_summary <- all_summary %>%
      mutate(feature_id = paste(chrom, start, end, sep = "_"))
    
    h2_cats <- c("Heritable", "Non-heritable", "Low prediction")
    
    for (h2_cat in h2_cats) {
      
      vmr_combined_venn <- all_summary %>%
        filter(h2_category == h2_cat) %>%
        group_by(r2_threshold) %>%
        summarise(regions = list(feature_id), .groups = "drop") %>%
        deframe()
      
      p <- ggvenn(vmr_combined_venn, fill_color = c("blue", "red"), 
                  fill_alpha = 0.4, stroke_size = 0, set_name_size = 5,
                  show_stats = "c", auto_scale = TRUE) +
        ggtitle(paste(tissue, h2_cat, pop, sep = " ")) +
        theme_void() + theme(legend.position = "none")
      
      plot_list[[paste(tissue, h2_cat, sep = "_")]] <- p
    }
  }
  combined_p <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 3)
  fn_venn <- file.path(out_path, paste0("combined_venn_", pop))
  save_plot(combined_p, fn_venn, 10, 10, 300)
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()