##### Performs Spearman's Correlation for Elastic-Net on Simulated Data. #####
suppressPackageStartupMessages({
    library(ggpubr)
    library(here)
    library(dplyr)
    library(readr)
    library(tidyr)
})

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
  }
}

out_path <- here("simulation-analysis/comparison/correlation/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

results_list <- list()
plot_list <- list()
num_indivs <- c(100, 150, 200, 250, 500, 1000, 5000)

for (num_indiv in num_indivs) {
  
  enet_path <- here("simulation-analysis/elastic-net", paste0("sim_", num_indiv, "_indiv"), "_m",
                    paste0("simulation_", num_indiv, "_summary_elastic-net.tsv"))
  target_path <- here("inputs", "simulated-data", "_m", paste0("sim_", num_indiv, "_indiv"), 
                      "snp_phenotype_mapping.tsv")
  
  enet <- read.table(enet_path, sep = "\t", header = TRUE)
  target <- read.table(target_path, sep = "\t", header = TRUE)
  
  merged <- inner_join(enet, target, by = c("pheno_id" = "phenotype_id"))

  merged <- merged %>%
    mutate(h2_category = case_when(
             r_squared_cv <= 0.75 ~ "Low prediction",
             target_heritability < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
             target_heritability > 0.1 & r_squared_cv > 0.75 ~ "Heritable"
           ))
  
  counts <- merged %>%
    group_by(h2_category) %>%
    summarise(n = n(), .groups = "drop")
  
  labels <- setNames(
    paste0(counts$h2_category, "\n(n = ", counts$n, ")"),
    counts$h2_category
  )
  
  if (nrow(merged) < 3) {
    warning(paste("Not enough data for N =", num_indiv, "- skipping"))
    next
  }
  
  # Spearman correlation
  spearman <- merged %>%
    group_by(h2_category) %>%
    summarise(
      spearman_rho = cor.test(target_heritability, h2_unscaled, method = "spearman")$estimate,
      p_value = cor.test(target_heritability, h2_unscaled, method = "spearman")$p.value,
      .groups = "drop"
    ) %>%
    mutate(N = num_indiv)
  
  results_list[[as.character(num_indiv)]] <- spearman
  
  # Plot for each sample size
  p <- ggscatter(merged, x = "target_heritability", y = "h2_unscaled",
                 add = "reg.line", size = 1, alpha = 0.5,
                 xlab = "True h²", ylab = "Estimated h²",
                 conf.int = TRUE,
                 cor.coef = TRUE, cor.coef.size = 4,
                 cor.method = "spearman",
                 color = "h2_category",
                 cor.coeff.args = list(
                   label.sep = "\n",
                   label.x.npc = 0.05,
                   label.y.npc = 0.95,
                   font.face = "bold"
                 ),
                 add.params = list(fill = "lightgray", alpha = 0.75),
                 ggtheme = theme_pubr(base_size = 15, border = TRUE)
  ) +
    facet_wrap(~h2_category, labeller = as_labeller(labels), scales = "free_x") +
    ggtitle(paste("N =", num_indiv)) +
    labs(color = NULL) +
    font("xy.title", face = "bold", size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "black")
  
  plot_list[[as.character(num_indiv)]] <- p
}

if (length(plot_list) < 8) {
  blank_plot <- ggplot() + theme_void()
  plot_list[["placeholder"]] <- blank_plot
}

# Combine all plots into a 2x3 grid
combined_plot <- ggarrange(plotlist = plot_list, ncol = 4, nrow = 2, labels = NULL)

# Save combined plot
plot_file <- file.path(out_path, "elastic_net_correlation_combined")
save_plot(combined_plot, plot_file, w = 20, h = 10)

# Combine and write results
results_df <- bind_rows(results_list)
print(results_df)

write.csv(results_df, file = file.path(out_path, "elastic_net_spearman_correlation_results.csv"), 
          row.names = FALSE)

# Reproducibility info
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
