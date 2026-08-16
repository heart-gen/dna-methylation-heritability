#### Plot distribution of heritability estimates ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
})

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

summarise_h2 <- function(vmr, tissue, pop, out_path) {
  summary_df <- vmr %>%
    group_by(h2_category) %>%
    summarise(
      n = n(),
      mean_h2 = mean(h2_unscaled, na.rm = TRUE),
      median_h2 = median(h2_unscaled, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_df)
  write.csv(summary_df, file = file.path(out_path, paste0("h2_summary_", 
                                         tissue, "_", pop, ".csv")), 
            row.names = FALSE)
}

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}

plot_density <- function(vmr, tissue) {
  counts <- vmr %>%
    group_by(h2_category) %>%
    summarise(n = n(), .groups = "drop")
  
  labels <- setNames(
    paste0(counts$h2_category, "\n(n = ", counts$n, ")"),
    counts$h2_category
  )

  category_colors <- c(
    "Heritable" = "#497C8A",
    "Non-heritable" = "#8CA77B",
    "Low prediction" = "#E3A27F"
  )
  
  tissue_title <- ifelse(tolower(tissue) == "dlpfc", "DLPFC", tools::toTitleCase(tissue))

  p_hist <- gghistogram(vmr, x = "h2_unscaled", 
                        add_density = TRUE, rug = TRUE, 
                        add = "mean",
                        color = "h2_category", fill = "h2_category",
                        ggtheme = theme_pubr(base_size = 20, border = TRUE),
                        xlab = "Estimated h²", ylab = "Count") +
    facet_wrap(~h2_category, labeller = as_labeller(labels), scales = "free_x") +
    scale_color_manual(values = category_colors) +
    scale_fill_manual(values = category_colors) +
    ggtitle(paste("h2 distribution:", tissue_title)) +
    labs(color = NULL, fill = NULL) +
    font("xy.title", face = "bold", size = 14) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "black") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p_hist)
}

# Main
tissues <- c("caudate", "hippocampus", "dlpfc")
pops <- c("AA_all", "EA_all", "AA_matched", "EA_matched")

out_path <- here("heritability/elastic_net_model/all_individuals/tissue_comparison/h2_distribution/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

for (pop in pops){

  hist_plots <- list()

  for (tissue in tissues) {

    if (pop %in% c("AA_all", "EA_all")){

      pop_fn <- gsub("_all", "", pop)

      # Read in summary table
      enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                        paste0(tissue, "/_m/", tissue, "_summary_elastic-net_", 
                        pop_fn, ".tsv"))
      enet <- read.table(enet_file, sep = "\t", header = TRUE)
      
      vmr <- filter_sites(enet)

    } else if (pop == "AA_matched") {

      # Read in summary table
      enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                        paste0(tissue, "/_m/", tissue, "_summary_elastic-net_matched_r2_0.3.tsv"))
      enet <- read.table(enet_file, sep = "\t", header = TRUE)

      vmr <- enet %>% na.omit() %>%
          select(chrom, start, end, tissue = region, 
                 r_squared_cv = r_squared_cv_AA,
                 h2_unscaled = h2_unscaled_AA, h2_category)

    } else if (pop == "EA_matched") {  

      # Read in summary table
      enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                        paste0(tissue, "/_m/", tissue, "_summary_elastic-net_matched_r2_0.3.tsv"))
      enet <- read.table(enet_file, sep = "\t", header = TRUE)

      vmr <- enet %>% na.omit() %>% 
          select(chrom, start, end, tissue = region, 
                 r_squared_cv = r_squared_cv_EA,
                 h2_unscaled = h2_unscaled_EA, h2_category)
    }

    # Summarize
    summarise_h2(vmr, tissue, pop, out_path)

    # Plotting
    p_hist <- plot_density(vmr, tissue)
      
    # Store plots
    hist_plots[[paste0(tissue, "_", pop)]] <- p_hist
  }

  # Save plots
  combined_hist <- ggarrange(plotlist = hist_plots, ncol = 3, nrow = 1)
  fn_hist <- file.path(out_path, paste0("vmr_h2_distribution_", pop))
  save_plot(combined_hist, fn_hist, 20, 6)

}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
