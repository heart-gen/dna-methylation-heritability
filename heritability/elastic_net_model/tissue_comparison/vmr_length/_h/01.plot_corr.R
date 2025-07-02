#### Correlate VMR length to heritability estimates ####

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
      r_squared_cv <= 0.75 ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
    ))
  return(vmr)
}

cal_vmr_length <- function(vmr) {
  vmr <- vmr %>%
    mutate(length = end - start)
  return(vmr)
}

get_long_vmrs <- function(vmr, tissue, out_path) {
  vmr_long <- vmr %>%
    filter(length > 5000) %>%
    select(chrom, start, end)
  
  write.table(vmr_long,
              file = file.path(out_path, paste0("long_vmrs_", tissue, ".bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  return(vmr_long)
}

summarise_length <- function(vmr, tissue, out_path) {
  summary_df <- vmr %>%
    group_by(h2_category) %>%
    summarise(
      n = n(),
      mean_length = mean(length, na.rm = TRUE),
      median_length = median(length, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_df)
  write.csv(summary_df, file = file.path(out_path, paste0("vmr_length_summary_", tissue, ".csv")), 
            row.names = FALSE)
}

spearman_corr <- function(vmr, tissue, out_path) {
  spearman <- vmr %>% 
    group_by(h2_category) %>%
    summarise(
      spearman_rho = cor.test(length, h2_unscaled, method = "spearman")$estimate,
      p_value = cor.test(length, h2_unscaled, method = "spearman")$p.value,
      n = n()
    )
  print(spearman)
  write.csv(spearman, file = file.path(out_path, paste0("vmr_length_h2_corr_", tissue, ".csv")), 
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
  
  p_hist <- gghistogram(vmr, x = "length", 
                        add_density = TRUE, rug = TRUE, 
                        add = "median",
                        color = "h2_category", fill = "h2_category",
                        ggtheme = theme_pubr(base_size = 15, border = TRUE),
                        xlab = "Length (BP)", ylab = "Count") +
    facet_wrap(~h2_category, labeller = as_labeller(labels), scales = "free_x") +
    ggtitle(paste("VMR distribution:", tissue)) +
    labs(color = NULL, fill = NULL) +
    font("xy.title", face = "bold", size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p_hist)
}

plot_corr <- function(vmr, tissue) {
  counts <- vmr %>%
    group_by(h2_category) %>%
    summarise(n = n(), .groups = "drop")
  
  labels <- setNames(
    paste0(counts$h2_category, "\n(n = ", counts$n, ")"),
    counts$h2_category
  )
  
  p_corr <- ggscatter(vmr, x = "length", y = "h2_unscaled",
                      add = "reg.line", size = 1, alpha = 0.5,
                      xlab = "VMR Length (BP)", ylab = "Estimated h2", 
                      conf.int = TRUE,
                      cor.coef = TRUE, cor.coef.size = 4,
                      cor.coeff.args = list(
                        label.sep = "\n",
                        label.x.npc = 0.05,
                        label.y.npc = 0.95,
                        font.face = "bold"
                      ),
                      cor.method = "spearman",
                      color = "h2_category",
                      add.params = list(fill = "lightgray", alpha = 0.75),
                      ggtheme = theme_pubr(base_size = 15, border = TRUE)
  ) +
    facet_wrap(~h2_category, labeller = as_labeller(labels), scales = "free_x") +
    labs(color = NULL) +
    ggtitle(paste("VMR length correlation:", tissue)) +
    font("xy.title", face = "bold", size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(p_corr)
}

# Main
tissues <- c("caudate", "hippocampus", "dlpfc")
out_path <- here("heritability/elastic_net_model/tissue_comparison/vmr_length/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

hist_plots <- list()
corr_plots <- list()

for (tissue in tissues) {
  # Read in summary table
  enet_file <- here("heritability/elastic_net_model/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
  enet <- read.table(enet_file, sep = "\t", header = TRUE)
  
  vmr <- filter_sites(enet)
  vmr <- cal_vmr_length(vmr)
  vmr_long <- get_long_vmrs(vmr, tissue, out_path)
  
  # Summarize
  summarise_length(vmr, tissue, out_path)

  # Spearman correlation test 
  spearman_corr(vmr, tissue, out_path)

  # Plotting
  p_hist <- plot_density(vmr, tissue)
  p_corr <- plot_corr(vmr, tissue)
  
  # Store plots
  hist_plots[[tissue]] <- p_hist
  corr_plots[[tissue]] <- p_corr
}

# Save plots
combined_hist <- ggarrange(plotlist = hist_plots, ncol = 3, nrow = 1)
combined_corr <- ggarrange(plotlist = corr_plots, ncol = 3, nrow = 1)

fn_hist <- file.path(out_path, "VMR_length_distribution")
fn_corr <- file.path(out_path, "VMR_length_h2_correlation")
save_plot(combined_hist, fn_hist, 20, 6)
save_plot(combined_corr, fn_corr, 20, 6)


#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
