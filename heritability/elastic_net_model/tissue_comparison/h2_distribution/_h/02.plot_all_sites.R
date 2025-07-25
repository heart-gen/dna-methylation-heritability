#### Correlate VMR length to heritability estimates ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
})

## Function
summarise_h2 <- function(vmr, out_path) {
  summary_df <- vmr %>%
    group_by(tissue) %>%
    summarise(
      n = n(),
      mean_h2 = mean(h2_unscaled, na.rm = TRUE),
      sd_h2 = sd(h2_unscaled, na.rm = TRUE),
      median_h2 = median(h2_unscaled, na.rm = TRUE),
      .groups = "drop"
    )
    write.csv(summary_df, file = file.path(out_path, paste0("h2_summary_all_sites.csv")), 
            row.names = FALSE)
  print(summary_df)
}

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}

plot_density <- function(vmr, tissue) {
  # Define palette
  tissue_colors <- c(
    "caudate" = "#7372A6",
    "dlpfc" = "#B36F61",
    "hippocampus" = "#C5AC47"
  )
  
  counts <- vmr %>%
    group_by(tissue) %>%
    summarise(n = n(), .groups = "drop")
  
  labels <- setNames(
    paste0(counts$tissue, "\n(n = ", counts$n, ")"),
    counts$tissue
  )
  
  p_hist <- gghistogram(vmr, x = "h2_unscaled", 
                        add_density = TRUE, rug = TRUE, 
                        add = "median",
                        color = "tissue", fill = "tissue",
                        ggtheme = theme_pubr(base_size = 15, border = TRUE),
                        xlab = "Estimated h²", ylab = "Count") +
    facet_wrap(~tissue, labeller = as_labeller(labels), scales = "free_x") +
    scale_color_manual(values = tissue_colors) +
    scale_fill_manual(values = tissue_colors) +
    ggtitle(paste("h2 distribution for all VMRs")) +
    labs(color = NULL, fill = NULL) +
    font("xy.title", face = "bold", size = 14) +
    geom_vline(xintercept = 0.19, linetype = "dashed", color = "black") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p_hist)
}

## Main
tissues <- c("caudate", "hippocampus", "dlpfc")

out_path <- here("heritability/elastic_net_model/tissue_comparison/h2_distribution/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

vmr_all    <- list()

for (tissue in tissues) {
  # Read in summary table
  enet_file <- here("heritability/elastic_net_model/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
  vmr <- read.table(enet_file, sep = "\t", header = TRUE)
    
  # Store vmrs across all tissues
  vmr$tissue <- tissue
  vmr_all[[tissue]] <- vmr
}

vmr_all <- bind_rows(vmr_all)

# Summarize h2 of vmrs
summarise_h2(vmr_all, out_path)

# Save plot
p_hist <- plot_density(vmr_all, tissue)
fn_hist <- file.path(out_path, "all_VMR_h2_distribution")
save_plot(p_hist, fn_hist, 14, 6)
