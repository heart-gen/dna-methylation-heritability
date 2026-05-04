#### Plot overlap of VMRs in the discovery and replication cohort ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggvenn)
  library(ggplot2)
  library(tidyverse)
  library(ggpubr)
})

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

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}
## Main
tissues <- c("caudate", "hippocampus", "dlpfc")

out_path <- here("heritability/elastic_net_model/all_individuals/tissue_comparison/venn_diagram/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

plot_list <- list()

for (tissue in tissues) {
  # Read in summary table
  enet_file_AA <- here("heritability/elastic_net_model/all_individuals/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net_AA.tsv"))
  enet_AA <- read.table(enet_file_AA, sep = "\t", header = TRUE)
  
  enet_file_EA <- here("heritability/elastic_net_model/all_individuals/", 
                       paste0(tissue, "/_m/", tissue, "_summary_elastic-net_EA.tsv"))
  enet_EA <- read.table(enet_file_EA, sep = "\t", header = TRUE)
  
  vmr_AA <- filter_sites(enet_AA)
  vmr_EA <- filter_sites(enet_EA)

  # Get matched VMRs by h2 cat
  vmr_all_indiv <- inner_join(vmr_AA, vmr_EA, by = c("chrom", "start", "end"),
                             suffix = c("_AA", "_EA")) %>%
    filter(h2_category_AA == h2_category_EA) %>%
    rename(h2_category = h2_category_AA) %>%
    select(-c(h2_category_EA, race_AA, race_EA, region_EA)) %>%
    mutate(feature_id = paste(chrom, start, end, sep = "_"),
           cohort = "replication")

  # Get BA only VMRs
  enet_file_BA_only <- here("heritability/elastic_net_model/BA_only/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
  enet_BA_only <- read.table(enet_file_BA_only, sep = "\t", header = TRUE)

  vmr_BA_only <- filter_sites(enet_BA_only) %>%
    mutate(feature_id = paste(chrom, start, end, sep = "_"),
           cohort = "discovery")

  vmr_all <- bind_rows(vmr_all_indiv, vmr_BA_only)
  
  h2_cats <- c("Heritable", "Non-heritable", "Low prediction")
  
  for (h2_cat in h2_cats) {
    
    vmr_combined_venn <- vmr_all %>%
      filter(h2_category == h2_cat) %>%
      group_by(cohort) %>%
      summarise(regions = list(feature_id), .groups = "drop") %>%
      deframe()
    
    p <- ggvenn(vmr_combined_venn, fill_color = c("brown", "blue"), 
                fill_alpha = 0.4, stroke_size = 0, set_name_size = 5,
                show_stats = "c", auto_scale = TRUE) +
      ggtitle(paste(tissue, h2_cat, sep = " ")) +
      theme_void() + theme(legend.position = "none")
    
    plot_list[[paste(tissue, h2_cat, sep = "_")]] <- p
  }
}

combined_p <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 3)
fn_venn <- file.path(out_path, "cohort_venn")
save_plot(combined_p, fn_venn, 10, 10)
  
#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()