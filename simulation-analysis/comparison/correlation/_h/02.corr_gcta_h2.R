##### Performs Spearman's Correlation for GCTA on Simulated Data. #####
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
num_indivs <- c(100, 150, 200, 250, 500, 1000, 5000, 10000)

for (num_indiv in num_indivs) {
  
  gcta_path <- here("simulation-analysis/gcta", "_m", "summary", "greml_summary.tsv")
  target_path <- here("inputs", "simulated-data", "_m", paste0("sim_", num_indiv, "_indiv"), 
                      "snp_phenotype_mapping.tsv")
  
  gcta <- read.table(gcta_path, sep = "\t", header = TRUE)
  target <- read.table(target_path, sep = "\t", header = TRUE)

  merged <- inner_join(gcta, target, by = c("pheno" = "phenotype_id"))
  filtered <- dplyr::filter(merged, N == num_indiv)

  filtered <- filtered %>%
    mutate(p_adjusted_fdr = p.adjust(Pval_Variance, method = "fdr"),
      h2_category = case_when(
          p_adjusted_fdr >= 0.05 ~ "Low prediction",
          target_heritability < 0.1 & p_adjusted_fdr < 0.05 ~ "Non-heritable",
          target_heritability > 0.1 & p_adjusted_fdr < 0.05 ~ "Heritable"
    ))
  
  counts <- filtered %>%
    group_by(h2_category) %>%
    summarise(n = n(), .groups = "drop")
  
  labels <- setNames(
    paste0(counts$h2_category, "\n(n = ", counts$n, ")"),
    counts$h2_category
  )
  
  if (nrow(filtered) < 3) {
    warning(paste("Not enough data for N =", num_indiv, "- skipping"))
    next
  }
  
  # Spearman correlation
  spearman <- filtered %>%
    group_by(h2_category) %>%
    summarise(
      spearman_rho = cor.test(target_heritability, Sum.of.V.G._Vp_Variance, method = "spearman")$estimate,
      p_value = cor.test(target_heritability, Sum.of.V.G._Vp_Variance, method = "spearman")$p.value,
      .groups = "drop"
    ) %>%
    mutate(N = num_indiv)
  
  results_list[[as.character(num_indiv)]] <- spearman
  
  # Plot for each sample size
  p <- ggscatter(filtered, x = "target_heritability", y = "Sum.of.V.G._Vp_Variance",
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  plot_list[[as.character(num_indiv)]] <- p
}
# Combine all plots into a 2x2 grid
combined_plot <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 2, labels = NULL)

# Save combined plot
plot_file <- file.path(out_path, "gcta_correlation_combined")
save_plot(combined_plot, plot_file, w = 18, h = 10)

# Combine and write results
results_df <- bind_rows(results_list)
print(results_df)

write.csv(results_df, file = file.path(out_path, "gcta_spearman_correlation_results.csv"), row.names = FALSE)

# Reproducibility info
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()