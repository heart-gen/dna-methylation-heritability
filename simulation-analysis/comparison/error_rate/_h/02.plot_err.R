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
         y = "Value",
         color = "Metric",
         shape = "Method",
         linetype = "Method") +
    theme_minimal(base_size = 15) +
    theme(legend.position = "right", 
          panel.grid.minor = element_blank())
  
  return(p)
}

## Main
# Create output directory
out_path <- here("simulation-analysis/comparison/error_rate/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# Read in error rates
error_rate_fn <- here("simulation-analysis/comparison/error_rate/_m/power-analysis.tsv")
error_rate <- read.table(error_rate_fn, sep = "\t", header = TRUE)
error_rate <- filter_error(error_rate)

# Plot
p <- plot_error(error_rate)

# Save plot
plot_file <- file.path(out_path, "simulated_error_rate")
save_plot(p, plot_file, w = 10, h = 10, dpi = 300)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
