#### Correlate VMR length to heritability estimates ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
})

## Function

cal_vmr_length <- function(vmr) {
  vmr <- vmr %>%
    mutate(length = end - start,
           log10_length = log10(length))
  return(vmr)
}

summarise_length <- function(vmr) {
  summary_df <- vmr %>%
    group_by(tissue) %>%
    summarise(
      n = n(),
      mean_length = mean(length, na.rm = TRUE),
      sd_length = sd(length, na.rm = TRUE),
      median_length = median(length, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_df)
}

spearman_corr <- function(vmr) {
  spearman <- vmr %>% 
    group_by(tissue) %>%
    summarise(
      spearman_rho = cor.test(log10_length, h2_unscaled, method = "spearman")$estimate,
      spearman_p_value = cor.test(log10_length, h2_unscaled, method = "spearman")$p.value,
      n = n()
    )
  print(spearman)
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
  
  p_hist <- gghistogram(vmr, x = "log10_length", 
                        add_density = TRUE, rug = TRUE, 
                        add = "median",
                        color = "tissue", fill = "tissue",
                        ggtheme = theme_pubr(base_size = 15, border = TRUE),
                        xlab = "Length (BP)", ylab = "Count") +
    facet_wrap(~tissue, labeller = as_labeller(labels), scales = "free_x") +
    scale_color_manual(values = tissue_colors) +
    scale_fill_manual(values = tissue_colors) +
    ggtitle(paste("Length distribution for all VMRs")) +
    labs(color = NULL, fill = NULL) +
    font("xy.title", face = "bold", size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p_hist)
}

plot_corr <- function(vmr, tissue) {
  counts <- vmr %>%
    group_by(tissue) %>%
    summarise(n = n(), .groups = "drop")
  
  labels <- setNames(
    paste0(counts$tissue, "\n(n = ", counts$n, ")"),
    counts$tissue
  )
  
  # Define palette
  tissue_colors <- c(
    "caudate" = "#7372A6",
    "dlpfc" = "#B36F61",
    "hippocampus" = "#C5AC47"
  )
  
  p_corr <- ggscatter(vmr, x = "log10_length", y = "h2_unscaled",
                      add = "reg.line", size = 1, alpha = 0.75,
                      xlab = "VMR Length (BP)", ylab = "Estimated h2", 
                      conf.int = TRUE,
                      cor.coef = TRUE, cor.coef.size = 4,
                      cor.coeff.args = list(
                        label.sep = "\n",
                        label.x.npc = 0.05,
                        label.y.npc = 0.95
                      ),
                      cor.method = "spearman",
                      color = "tissue",
                      add.params = list(fill = "lightgray", alpha = 0.75),
                      ggtheme = theme_pubr(base_size = 15, border = TRUE)
  ) +
    facet_wrap(~tissue, labeller = as_labeller(labels), scales = "free_x") +
    scale_color_manual(values = tissue_colors) +
    scale_fill_manual(values = tissue_colors) +
    labs(color = NULL) +
    ggtitle(paste("Correlation of H2 and VMR length for all sites")) +
    font("xy.title", face = "bold", size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "#2A0F07")
  return(p_corr)
}

## Main
tissues <- c("caudate", "hippocampus", "dlpfc")

out_path <- here("heritability/elastic_net_model/tissue_comparison/vmr_length/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

vmr_all    <- list()

for (tissue in tissues) {
  # Read in summary table
  enet_file <- here("heritability/elastic_net_model/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
  vmr <- read.table(enet_file, sep = "\t", header = TRUE)
  
  vmr <- cal_vmr_length(vmr)
  
  # Store vmrs across all tissues
  vmr$tissue <- tissue
  vmr_all[[tissue]] <- vmr
}

vmr_all <- bind_rows(vmr_all)

# Summarize length of vmrs
summarise_length(vmr_all)

# Spearman correlation test 
spearman_corr(vmr_all)

p_hist <- plot_density(vmr_all, tissue)
p_corr <- plot_corr(vmr_all, tissue)

fn_hist <- file.path(out_path, "all_VMR_length_distribution")
fn_corr <- file.path(out_path, "all_VMR_length_h2_correlation")
save_plot(p_hist, fn_hist, 14, 6)
save_plot(p_corr, fn_corr, 14, 6)