#### Plot error rate ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
  library(tidyverse)
})

## Function 
save_plot <- function(p, fn, w, h, dpi){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h, dpi=dpi)
  }
}

filter_error <- function(error_rate) {
  # Remove failed sample sizes
  filtered_error_rate <- error_rate %>%
    filter(!(method == "gcta" & sample_size %in% c(100, 150, 200, 250)))
  
  # Convert to long format for plotting
  error_rate <- filtered_error_rate %>%
    pivot_longer(cols = c(power, type1_error, type2_error),
                 names_to = "metric",
                 values_to = "value")
  
  error_rate <- error_rate %>%
    mutate(
      method = fct_recode(method,
                          "GREML-LDMS" = "gcta",
                          "Elastic-net" = "elastic-net"),
      metric = fct_recode(metric,
                          "Type 1 Error" = "type1_error",
                          "Type 2 Error" = "type2_error",
                          "Power" = "power")
    )
  return(error_rate)
}

plot_error <- function(error_rate) {
  p <- ggplot(error_rate, aes(x = sample_size, y = value, color = metric, shape = method)) +
    geom_point(size = 5) +
    geom_line(aes(linetype = method), size = 1) +
    scale_x_continuous(trans = "log10", breaks = unique(error_rate$sample_size)) +
    scale_color_manual(
      values = c("Power" = "#7B8C99",
                 "Type 1 Error" = "#B35A4E",
                 "Type 2 Error" = "#D4BFAA")
    ) +
    labs(x = "Sample Size (N)",
         y = "Error Measurement",
         color = "Metric",
         shape = "Method",
         linetype = "Method") +
    theme_minimal(base_size = 20) +
    theme(legend.position = "right", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 1, linetype = "solid"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

plot_bar <- function(error_rate_all, n_samples) {
  error_rate <- error_rate_all %>% 
    filter(sample_size == n_samples, method == "Elastic-net") %>%
    mutate(r2_threshold = as.factor(r2_threshold))
  
  p_bar <- ggplot(error_rate, aes(x = r2_threshold, y = value, fill = metric)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(
      values = c("Power" = "#7B8C99",
                 "Type 1 Error" = "#B35A4E",
                 "Type 2 Error" = "#D4BFAA")
    ) +
    labs(x = "r2 Threshold",
         y = "Error Measurement",
         fill = "Metric",
         title = "Error Estimates at N = 250") +
    theme_minimal(base_size = 20) +
    theme(legend.position = "right", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 1, linetype = "solid"),
          axis.text.x = element_text(angle = 45, hjust = 1))

  return(p_bar)
}

## Main
# Create output directory
out_path <- here("simulation-analysis/comparison/error_rate/r2_testing/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# Get r2 threshold from command line
args <- commandArgs(trailingOnly = TRUE)
r2_thresh <- as.numeric(args[1])

# Read in error rates
error_rate_fn <- here("simulation-analysis/comparison/error_rate/r2_testing/_m", paste0("power-analysis_r2_", r2_thresh, ".tsv"))
error_rate <- read.table(error_rate_fn, sep = "\t", header = TRUE)
error_rate <- filter_error(error_rate)

# Combine error rates for plotting
r2_vals <- c(0.3, 0.5, 0.7)

error_list <- lapply(r2_vals, function(r2_thresh) {
  
  error_rate_fn <- here("simulation-analysis/comparison/error_rate/r2_testing/_m",
             paste0("power-analysis_r2_", r2_thresh, ".tsv"))
  
  error_rate <- read.table(error_rate_fn, sep = "\t", header = TRUE)
  error_rate$r2_threshold <- r2_thresh
  
  filter_error(error_rate)
})

error_rate_all <- bind_rows(error_list)

# Plot
p <- plot_error(error_rate)
n_samples <- 250
p_bar <- plot_bar(error_rate_all, n_samples)

# Save plots
plot_file <- file.path(out_path, paste0("simulated_error_rate_r2_", r2_thresh))
barplot_file <- file.path(out_path, paste0("simulated_error_rate_n_", n_samples))
save_plot(p_bar, barplot_file, w = 6, h = 4, dpi = 300)
save_plot(p, plot_file, w = 10, h = 6, dpi = 300)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
