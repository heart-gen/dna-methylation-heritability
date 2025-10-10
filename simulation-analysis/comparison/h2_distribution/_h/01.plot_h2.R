#### Plot distribution of heritability estimates ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
})

filter_sites <- function(enet) {
  enet <- na.omit(enet)
  
  enet_all <- enet %>% 
    mutate(h2_category = "All sites")

  enet_heritable <- enet %>% 
    filter(h2_unscaled > 0.1, r_squared_cv >= 0.75) %>%
    mutate(h2_category = "Heritable")

  enet_non_heritable <- enet %>% 
    filter(h2_unscaled < 0.1, r_squared_cv >= 0.75) %>%
    mutate(h2_category = "Non-Heritable")

  enet_low_pred <- enet %>% 
    filter(r_squared_cv < 0.75) %>%
    mutate(h2_category = "Low Prediction")

  # Combine all for plotting
  enet_combined <- bind_rows(
    enet_all,
    enet_heritable,
    enet_non_heritable,
    enet_low_pred
  )
  return(enet_combined)
}

summarise_h2 <- function(enet_combined, num_indiv, out_path) {
  summary_df <- enet_combined %>%
    group_by(h2_category) %>%
    summarise(
      n = n(),
      mean_h2 = mean(h2_unscaled, na.rm = TRUE),
      median_h2 = median(h2_unscaled, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_df)
  write.csv(summary_df, file = file.path(out_path, paste0("h2_summary_", num_indiv, ".csv")), 
            row.names = FALSE)
}

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}

plot_density <- function(enet_combined, num_indiv) {
  counts <- enet_combined %>%
    group_by(h2_category) %>%
    summarise(n = n(), .groups = "drop")
  
  labels <- setNames(
    paste0(counts$h2_category, "\n(n = ", counts$n, ")"),
    counts$h2_category
  )
  
  category_colors <- c(
    "All sites" = "#7B8C99",
    "Heritable" = "#497C8A",
    "Non-heritable" = "#8CA77B",
    "Low prediction" = "#E3A27F"
  )

  p_hist <- gghistogram(enet_combined, x = "h2_unscaled", 
                        add_density = TRUE, rug = TRUE, 
                        add = "median",
                        color = "h2_category", fill = "h2_category",
                        ggtheme = theme_pubr(base_size = 15),
                        xlab = "Estimated h²", ylab = "Count") +
    ggtitle(paste("N =", num_indiv)) +
    labs(color = "h2 Category", fill = "h2 Category") +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "black") +
    scale_fill_manual(values = category_colors) +
    scale_color_manual(values = category_colors) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  return(p_hist)
}

# Main
out_path <- here("simulation-analysis/comparison/h2_distribution/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

hist_plots <- list()
num_indivs <- c(100, 150, 200, 250, 500, 1000, 5000)

for (num_indiv in num_indivs) {
  # Read in summary table
  enet_file <- here("simulation-analysis/elastic-net", paste0("sim_", num_indiv, "_indiv"), "_m",
                    paste0("simulation_", num_indiv, "_summary_elastic-net.tsv"))
  enet <- read.table(enet_file, sep = "\t", header = TRUE)
  
  enet_combined <- filter_sites(enet)
  
  # Summarize
  summarise_h2(enet_combined, as.character(num_indiv), out_path)
  
  # Plotting
  p_hist <- plot_density(enet_combined, as.character(num_indiv))
  
  # Store plots
  hist_plots[[as.character(num_indiv)]] <- p_hist
}

combined_hist <- ggarrange(plotlist = hist_plots, ncol = 4, nrow = 2)
fn_hist <- file.path(out_path, "enet_h2_distribution_overlap")
save_plot(combined_hist, fn_hist, 20, 10)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
