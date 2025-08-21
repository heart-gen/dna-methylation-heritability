#### Correlate number of SNPs to VMR length and heritability estimates ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
  library(rstatix)
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

cal_vmr_length <- function(vmr) {
  vmr <- vmr %>%
    mutate(length = end - start,
           log10_length = log10(length))
  return(vmr)
}

summarise_num_snps <- function(vmr, tissue, out_path) {
  summary_df <- vmr %>%
    group_by(h2_category) %>%
    summarise(
      n = n(),
      mean_num_snps = mean(num_snps, na.rm = TRUE),
      median_num_snps = median(num_snps, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_df)
  write.csv(summary_df, 
            file = file.path(out_path, 
                             paste0("num_snps_summary_", tissue, ".csv")), 
            row.names = FALSE)
}

spearman_corr <- function(vmr, tissue, out_path) {
  spearman <- vmr %>% 
    group_by(h2_category) %>%
    summarise(
      cor_snps_length      = cor.test(num_snps, length, method = "spearman")$estimate,
      p_snps_length        = cor.test(num_snps, length, method = "spearman")$p.value,
      cor_snps_h2          = cor.test(num_snps, h2_unscaled, method = "spearman")$estimate,
      p_snps_h2            = cor.test(num_snps, h2_unscaled, method = "spearman")$p.value,
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      p_snps_length_adj = p.adjust(p_snps_length, method = "fdr"),
      p_snps_h2_adj     = p.adjust(p_snps_h2, method = "fdr")
    )
  print(spearman)
  write.csv(spearman, 
            file = file.path(out_path, 
                             paste0("num_snps_corr_", tissue, ".csv")), 
            row.names = FALSE)
}

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}

plot_density <- function(vmr, tissue) {
  # Define palette
  category_colors <- c(
    "Heritable" = "#497C8A",
    "Non-heritable" = "#8CA77B",
    "Low prediction" = "#E3A27F"
  )
  
  counts <- vmr %>%
    group_by(h2_category) %>%
    summarise(n = n(), .groups = "drop")
  
  labels <- setNames(
    paste0(counts$h2_category, "\n(n = ", counts$n, ")"),
    counts$h2_category
  )
  
  tissue_title <- ifelse(tolower(tissue) == "dlpfc", "DLPFC", tools::toTitleCase(tissue))
  
  p_hist <- gghistogram(vmr, x = "num_snps", 
                        add_density = TRUE, rug = TRUE, 
                        add = "median",
                        color = "h2_category", fill = "h2_category",
                        ggtheme = theme_pubr(base_size = 15, border = TRUE),
                        xlab = "Number of SNPs", ylab = "Count") +
    facet_wrap(~h2_category, labeller = as_labeller(labels), scales = "free_x") +
    scale_color_manual(values = category_colors) +
    scale_fill_manual(values = category_colors) +
    ggtitle(paste("Number of SNPs distribution:", tissue_title)) +
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
    group_by(h2_category) %>%
    summarise(n = n(), .groups = "drop")
  
  labels <- setNames(
    paste0(counts$h2_category, "\n(n = ", counts$n, ")"),
    counts$h2_category
  )
  
  # Define palette
  category_colors <- c(
    "Heritable" = "#497C8A",
    "Non-heritable" = "#8CA77B",
    "Low prediction" = "#E3A27F"
  )
  
  tissue_title <- ifelse(tolower(tissue) == "dlpfc", "DLPFC", tools::toTitleCase(tissue))
  
  p_corr <- ggscatter(vmr, x = "num_snps", y = "h2_unscaled",
                      add = "reg.line", size = 1, alpha = 0.75,
                      xlab = "Number of SNPs", ylab = "Estimated h2", 
                      conf.int = TRUE,
                      cor.coef = TRUE, cor.coef.size = 4,
                      cor.coeff.args = list(
                        label.sep = "\n",
                        label.x.npc = 0.05,
                        label.y.npc = 0.95
                      ),
                      cor.method = "spearman",
                      color = "h2_category",
                      add.params = list(fill = "lightgray", alpha = 0.75),
                      ggtheme = theme_pubr(base_size = 15, border = TRUE)
  ) +
    facet_wrap(~h2_category, labeller = as_labeller(labels), scales = "free_x") +
    scale_color_manual(values = category_colors) +
    scale_fill_manual(values = category_colors) +
    labs(color = NULL) +
    ggtitle(paste("Number of SNPs vs. Estimated H2:", tissue_title)) +
    font("xy.title", face = "bold", size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "#2A0F07")
  return(p_corr)
}

plot_corr_snps <- function(vmr, tissue) {
  counts <- vmr %>%
    group_by(h2_category) %>%
    summarise(n = n(), .groups = "drop")
  
  labels <- setNames(
    paste0(counts$h2_category, "\n(n = ", counts$n, ")"),
    counts$h2_category
  )
  
  # Define palette
  category_colors <- c(
    "Heritable" = "#497C8A",
    "Non-heritable" = "#8CA77B",
    "Low prediction" = "#E3A27F"
  )
  
  tissue_title <- ifelse(tolower(tissue) == "dlpfc", "DLPFC", tools::toTitleCase(tissue))
  
  p_corr_snp <- ggscatter(vmr, x = "num_snps", y = "log10_length",
                      add = "reg.line", size = 1, alpha = 0.75,
                      xlab = "Number of SNPs", ylab = "log10(VMR Length)", 
                      conf.int = TRUE,
                      cor.coef = TRUE, cor.coef.size = 4,
                      cor.coeff.args = list(
                        label.sep = "\n",
                        label.x.npc = 0.05,
                        label.y.npc = 0.95
                      ),
                      cor.method = "spearman",
                      color = "h2_category",
                      add.params = list(fill = "lightgray", alpha = 0.75),
                      ggtheme = theme_pubr(base_size = 15, border = TRUE)
  ) +
    facet_wrap(~h2_category, labeller = as_labeller(labels), scales = "free_x") +
    scale_color_manual(values = category_colors) +
    scale_fill_manual(values = category_colors) +
    labs(color = NULL) +
    ggtitle(paste("Number of SNPs vs. VMR length:", tissue_title)) +
    font("xy.title", face = "bold", size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  return(p_corr_snp)
}

plot_violin <- function(vmr_all, snp_test) {
  heritability_colors <- c(
    "Heritable" = "#497C8A",
    "Non-heritable" = "#8CA77B",
    "Low prediction" = "#E3A27F"
  )
  
  p_violin <- ggviolin(vmr_all, x = "h2_category",
                             y = "num_snps", fill = "h2_category", 
                             color = "h2_category", facet.by = "tissue_title", 
                             add = c("jitter", "median"), 
                             alpha = 0.5, palette = heritability_colors, 
                             trim = FALSE
  ) +
    stat_pvalue_manual(snp_test, label = "p.adj.signif", tip.length = 0.01) +
    theme_pubr(base_size = 15, border = TRUE) +
    labs(
      title = "Number of SNP Differences Across Heritability Categories",
      x = "Heritability Category", y = "Number of SNPs"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  return(p_violin)
  
}

## Main
tissues <- c("caudate", "hippocampus", "dlpfc")

out_path <- here("heritability/elastic_net_model/tissue_comparison/vmr_length/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

hist_plots     <- list()
corr_plots     <- list()
corr_plots_snp <- list()

for (tissue in tissues) {
  # Read in summary table
  enet_file <- here("heritability/elastic_net_model/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
  enet <- read.table(enet_file, sep = "\t", header = TRUE)
  
  vmr <- filter_sites(enet)
  vmr <- cal_vmr_length(vmr)
  vmr$tissue <- tissue
  
  vmr_all <- bind_rows(if (exists("vmr_all")) vmr_all else NULL, vmr)
  vmr_all <- vmr_all %>%
    mutate(tissue_title = ifelse(tolower(tissue) == "dlpfc", "DLPFC", tools::toTitleCase(tissue)))
  
  # Summarize length of vmrs
  summarise_num_snps(vmr, tissue, out_path)
  
  # Spearman correlation test 
  spearman_corr(vmr, tissue, out_path)
  
  # Plotting
  p_hist <- plot_density(vmr, tissue)
  p_corr <- plot_corr(vmr, tissue)
  p_corr_snp <- plot_corr_snps(vmr, tissue)

  # Store plots
  hist_plots[[tissue]] <- p_hist
  corr_plots[[tissue]] <- p_corr
  corr_plots_snp[[tissue]] <- p_corr_snp
}

snp_test <- vmr_all %>%
  group_by(tissue_title) %>%
  pairwise_wilcox_test(
    num_snps ~ h2_category,
    p.adjust.method = "fdr"
  ) %>%
  add_y_position()

p_violin <- plot_violin(vmr_all, snp_test)

# Save plots
combined_hist <- ggarrange(plotlist = hist_plots, ncol = 3, nrow = 1)
combined_corr <- ggarrange(plotlist = corr_plots, ncol = 3, nrow = 1)
combined_corr_snp <- ggarrange(plotlist = corr_plots_snp, ncol = 3, nrow = 1)

fn_hist <- file.path(out_path, "num_snp_distribution")
fn_corr <- file.path(out_path, "num_snp_h2_correlation")
fn_corr_snp <- file.path(out_path, "num_snp_VMR_length_correlation")
fn_viol <- file.path(out_path, "num_snp_sig")

save_plot(combined_hist, fn_hist, 20, 6)
save_plot(combined_corr, fn_corr, 20, 6)
save_plot(combined_corr_snp, fn_corr_snp, 20, 6)
save_plot(p_violin, fn_viol, 14, 8)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()