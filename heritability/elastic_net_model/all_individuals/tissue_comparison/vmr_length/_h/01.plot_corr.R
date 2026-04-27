#### Correlate VMR length to heritability estimates ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
})

## Function
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

cal_vmr_length <- function(vmr) {
  vmr <- vmr %>%
    mutate(length = end - start,
           log10_length = log10(length))
  return(vmr)
}

get_long_vmrs <- function(vmr, tissue, pop, out_path) {
  vmr_long <- vmr %>%
    filter(length > 5000) %>%
    select(chrom, start, end)
  
  write.table(vmr_long,
              file = file.path(out_path, paste0("long_vmrs_", tissue, 
                               "_", pop, ".bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  return(vmr_long)
}

summarise_length <- function(vmr, tissue, pop, out_path) {
  summary_df <- vmr %>%
    group_by(h2_category) %>%
    summarise(
      n = n(),
      mean_length = mean(length, na.rm = TRUE),
      median_length = median(length, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_df)
  write.csv(summary_df, 
            file = file.path(out_path, 
                             paste0("vmr_length_summary_", tissue, "_", 
                                    pop, ".csv")), 
            row.names = FALSE)
}

spearman_corr <- function(vmr, tissue, pop, out_path) {
  spearman <- vmr %>% 
    group_by(h2_category) %>%
    summarise(
      spearman_rho = cor.test(log10_length, h2_unscaled, method = "spearman")$estimate,
      spearman_p_value = cor.test(log10_length, h2_unscaled, method = "spearman")$p.value,
      n = n()
    )
  print(spearman)
  write.csv(spearman, 
            file = file.path(out_path, 
                             paste0("vmr_length_h2_corr_", tissue, "_", 
                                    pop, ".csv")), 
            row.names = FALSE)
}

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}

plot_density <- function(vmr, tissue, pop) {
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
  
  p_hist <- gghistogram(vmr, x = "log10_length", 
                        add_density = TRUE, rug = TRUE, 
                        add = "median",
                        color = "h2_category", fill = "h2_category",
                        ggtheme = theme_pubr(base_size = 15, border = TRUE),
                        xlab = "log10(VMR Length)", ylab = "Count") +
    facet_wrap(~h2_category, labeller = as_labeller(labels), scales = "free_x") +
    scale_color_manual(values = category_colors) +
    scale_fill_manual(values = category_colors) +
    ggtitle(paste0("VMR length distribution:", tissue_title, "-", pop)) +
    labs(color = NULL, fill = NULL) +
    font("xy.title", face = "bold", size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p_hist)
}

plot_corr <- function(vmr, tissue, pop) {
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
  
  p_corr <- ggscatter(vmr, x = "log10_length", y = "h2_unscaled",
                      add = "reg.line", size = 1, alpha = 0.75,
                      xlab = "log10(VMR Length)", ylab = "Estimated h2", 
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
    ggtitle(paste0("VMR length correlation:", tissue_title, "-", pop)) +
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
pops <- c("AA", "EA")

out_path <- here("heritability/elastic_net_model/all_individuals/tissue_comparison/vmr_length/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

hist_plots <- list()
corr_plots <- list()

for (tissue in tissues) {
  for (pop in pops){

    # Read in summary table
    enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net_", 
                    pop, ".tsv"))
    enet <- read.table(enet_file, sep = "\t", header = TRUE)

    vmr <- filter_sites(enet)
    vmr <- cal_vmr_length(vmr)

    # Save VMRs > 5000bp
    vmr_long <- get_long_vmrs(vmr, tissue, pop, out_path)
  
    # Summarize length of vmrs
    summarise_length(vmr, tissue, pop, out_path)

    # Spearman correlation test 
    spearman_corr(vmr, tissue, pop, out_path)

    # Plotting
    p_hist <- plot_density(vmr, tissue, pop)
    p_corr <- plot_corr(vmr, tissue, pop)
    
    # Store plots
    hist_plots[[paste(pop, tissue)]] <- p_hist
    corr_plots[[paste(pop, tissue)]] <- p_corr
  }
}

# Save plots
combined_hist <- ggarrange(plotlist = hist_plots, ncol = 3, nrow = 2)
combined_corr <- ggarrange(plotlist = corr_plots, ncol = 3, nrow = 2)

fn_hist <- file.path(out_path, "VMR_length_distribution")
fn_corr <- file.path(out_path, "VMR_length_h2_correlation")
save_plot(combined_hist, fn_hist, 18, 14)
save_plot(combined_corr, fn_corr, 18, 14)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
