# Summarize simulated data using stacked barplot

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
  library(tidyr)
})

## Function
save_plot <- function(p, fn, w, h, dpi){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn,ext), plot=p, width=w, height=h, dpi=dpi)
  }
}

## Main
# Create output directory
out_path <- here("simulation-analysis/comparison/h2_distribution/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# Read in and filter summary tables
num_indivs <- c(100, 250, 500, 5000)
enet_list <- list()
target_list <- list()

for (num_indiv in num_indivs) {
  enet_file <- here("simulation-analysis/elastic-net", paste0("sim_", num_indiv, "_indiv"), "_m",
                    paste0("simulation_", num_indiv, "_summary_elastic-net.tsv"))
  enet <- read.table(enet_file, sep = "\t", header = TRUE)
  enet <- na.omit(enet)
  enet$N <- num_indiv
  enet_list[[as.character(num_indiv)]] <- enet
  
  target_path <- here("inputs", "simulated-data", "_m", paste0("sim_", num_indiv, "_indiv"), 
                      "snp_phenotype_mapping.tsv")
  target <- read.table(target_path, sep = "\t", header = TRUE)
  target$N <- num_indiv
  target_list[[as.character(num_indiv)]] <- target
}

enet_all   <- bind_rows(enet_list)
target_all <- bind_rows(target_list)

enet_all <- enet_all %>%
  mutate(method = "Elastic Net",
         h2_category = case_when(
           r_squared_cv <= 0.75 ~ "Low prediction",
           h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
           h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
         ))

target_all <- target_all %>%
  mutate(method = "Target",
         h2_category = case_when(
           target_heritability < 0.1 ~ "Non-heritable",
           target_heritability >= 0.1 ~ "Heritable"
         ))

gcta_path <- here("simulation-analysis/gcta", "_m", "summary", "greml_summary.tsv")
gcta <- read.table(gcta_path, sep = "\t", header = TRUE)
gcta <- gcta %>%
  mutate(method = "GREML-LDMS",
         p_adjusted_fdr = p.adjust(Pval_Variance, method = "fdr"),
         h2_category = case_when(
           p_adjusted_fdr >= 0.05 ~ "Low prediction",
           Sum.of.V.G._Vp_Variance < 0.1 & p_adjusted_fdr < 0.05 ~ "Non-heritable",
           Sum.of.V.G._Vp_Variance > 0.1 & p_adjusted_fdr < 0.05 ~ "Heritable"
         ))


enet_long <- enet_all %>%
  select(N, h2_unscaled, method) %>%
  rename(h2 = h2_unscaled)

target_long <- target_all %>%
  select(N, target_heritability, method) %>%
  rename(h2 = target_heritability)

gcta_long <- gcta %>%
  select(N, Sum.of.V.G._Vp_Variance, method) %>%
  rename(h2 = Sum.of.V.G._Vp_Variance)

# Combine all into one data frame
h2_combined <- bind_rows(target_long, enet_long, gcta_long)
h2_combined$method <- factor(h2_combined$method, levels = c("Target", "Elastic Net", "GREML-LDMS"))

method_colors <- c(
  "Target" = "#7B8C99",
  "Elastic Net" = "#B35A4E",
  "GREML-LDMS" = "#D4BFAA"
)

plot_list <- list()

for (n in num_indivs) {
  
  medians <- h2_combined %>%
    filter(N == n) %>%
    group_by(method) %>%
    summarize(median_h2 = median(h2, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(h2_combined %>% filter(N == n), aes(x = h2, fill = method, color = method)) +
    geom_density(alpha = 0.4, size = 1) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    labs(
      title = (paste("N =", n)),
      x = NULL,
      y = NULL,
      fill = "Method",
      color = "Method"
    ) +
    theme_minimal(base_size = 20) +
    geom_vline(xintercept = 0.1, linetype = "dashed", 
               color = "black", size = 1) +
    geom_vline(data = medians,
               aes(xintercept = median_h2, color = method),
               linetype = "dashed", size = 1) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
      
  plot_list[[as.character(n)]] <- p
  
}

combined_hist <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 2,
                           common.legend = TRUE, legend = "top")

combined_hist_annotated <- annotate_figure(
  combined_hist,
  left = text_grob("Count", rot = 90, size = 20),
  bottom = text_grob("Estimated heritability (h²)", size = 20)
)

fn_hist <- file.path(out_path, "simulated_h2_distribution")
save_plot(combined_hist_annotated, fn_hist, 8, 6, dpi = 300)

