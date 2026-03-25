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

get_summary <- function(enet_all, r2_threshold){
  enet_cat <- enet_all %>%
    mutate(h2_category = case_when(
             r_squared_cv <= r2_threshold ~ "Low prediction",
             h2_unscaled < 0.1 & r_squared_cv > r2_threshold ~ "Non-heritable",
             h2_unscaled >= 0.1 & r_squared_cv > r2_threshold ~ "Heritable"
           ),
           r2_threshold = r2_threshold)

  enet_cat <- enet_cat %>%
    mutate(h2_category = factor(h2_category,
                                levels = c("Heritable", "Non-heritable", "Low prediction")),
           region = recode(region, "caudate" = "Caudate"))

  return(enet_summary)
}

plot_venn <- function(enet_summary, tissue, h2_cat){
  
                                       # Prepare data for ggvenn
  vmr_combined_venn <- enet_summary %>%
    filter(h2_category == h2_cat) %>%
    group_by(r2_threshold) %>%
    summarise(regions = list(feature_id), .groups = "drop") %>%
    deframe()

                                        # Plot venn
  p <- ggvenn(vmr_combined_venn, fill_color = c("blue", "red"), 
              fill_alpha = 0.4, stroke_size = 0, set_name_size = 5,
              show_stats = "c", auto_scale = TRUE) +
    ggtitle(paste(tissue, h2_cat, sep = " ")) +
    theme_void() + theme(legend.position = "none")

  return(p)
}

## Main 
                                        # Create output directory 
out_path <- here("sensitivity-analysis/BA_only/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}
                                        # Define tissues
tissues <- c("caudate", "hippocampus", "dlpfc")
plot_list <- list()

for (tissue in tissues) {
                                        # Read in elastic net summaries
  enet_file <- here("heritability/elastic_net_model/BA_only/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
  enet <- read.table(enet_file, sep = "\t", header = TRUE)
  
                                        # Get summary for each r2 threshold
  enet_high <- get_summary(enet_all, r2_threshold = 0.75)
  enet_low <- get_summary(enet_all, r2_threshold = 0.3)

                                        # Combine summaries
  all_summary <- bind_rows(enet_high, enet_low) %>%
    mutate(feature_id = paste(chrom, start, end, sep = "_")) # Add for plotting
  
                                        # Define heritability categories
  h2_cats <- c("Heritable", "Non-heritable", "Low prediction")
  
  for (h2_cat in h2_cats) {
                                        # Plot venn for each tissue
                                        # and h2 category
    p <- plot_venn(all_summary, tissue, h2_cat)
    plot_list[[paste(tissue, h2_cat, sep = "_")]] <- p
  }
}

                                        # Create 3 x 3 grid 
combined_p <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 3)

                                        # Save plot
fn_venn <- file.path(out_path, "combined_venn")
save_plot(combined_p, fn_venn, 10, 10, 300)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
