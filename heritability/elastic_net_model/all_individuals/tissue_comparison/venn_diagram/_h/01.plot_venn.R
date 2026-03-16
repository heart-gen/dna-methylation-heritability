#### Plot overlap of VMRs in AA and EA ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggVennDiagram)
  library(ggplot2)
  library(tidyverse)
})

## Function
filter_sites <- function(enet) {
  vmr <- na.omit(enet)
  vmr <- vmr %>%
    mutate(h2_category = case_when(
      r_squared_cv <= 0.75 ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
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
  
  vmr_combined <- bind_rows(vmr_AA, vmr_EA)
    
  enet_summary <- vmr_combined %>%
    group_by(race, h2_category) %>%
    summarise(count = n(), .groups = "drop")
  print(enet_summary)
  
  vmr_combined <- vmr_combined %>%
    mutate(feature_id = paste(chrom, start, end, sep = "_"))
  
  h2_cats <- c("Heritable", "Non-heritable", "Low prediction")
  
  for (h2_cat in h2_cats) {
    
    vmr_combined_venn <- vmr_combined %>%
      filter(h2_category == h2_cat) %>%
      group_by(race) %>%
      summarise(regions = list(feature_id), .groups = "drop") %>%
      deframe()
    
    p <- ggVennDiagram(vmr_combined_venn) +
      ggtitle(paste0(tissue, ": ", h2_cat, " VMRs")) +
      theme(legend.position = "none")
    
    fn_venn <- file.path(out_path, paste0(tissue, "_", h2_cat, "_venn"))
    save_plot(p, fn_venn, 4, 4)
    
  }

}
  