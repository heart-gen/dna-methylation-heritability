##### Performs Spearman's Correlation for Elastic-Net on Simulated Data. #####
suppressPackageStartupMessages({
    library(ggplot2)
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

out_path <- here("simulation-analysis/elastic-net/correlation")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

results_list <- list()
num_indivs <- c(100, 150, 200, 250, 500, 1000)

for (num_indiv in num_indivs) {
  
  enet_path <- here("simulation-analysis/elastic-net", paste0("sim_", num_indiv, "_indiv"), "_m",
                    paste0("simulation_", num_indiv, "_summary_elastic-net.tsv"))
  target_path <- here("inputs", "simulated-data", "_m", paste0("sim_", num_indiv, "_indiv"), 
                      "snp_phenotype_mapping.tsv")
  
  enet <- read.table(enet_path, sep = "\t", header = TRUE)
  target <- read.table(target_path, sep = "\t", header = TRUE)

  merged <- inner_join(enet, target, by = c("pheno_id" = "phenotype_id"))
  merged <- dplyr::filter(merged, target_heritability > 0.1)
  
  # Spearman correlation
  rho <- cor(merged$h2_unscaled, merged$target_heritability, method = "spearman")
  
  results_list[[as.character(num_indiv)]] <- data.frame(
    N = num_indiv,
    spearman_rho = rho
  )
  
  # Plot for each sample size
  p <- ggplot(merged, aes(x = h2_unscaled, y = target_heritability)) +
    geom_point() +
    geom_smooth(method = "lm", color = "blue", se = FALSE) +
    labs(
      x = "Estimated h2",
      y = "Target h2",
      title = paste0("Simulated Correlation Results (Heritable sites, N = ", num_indiv, ")")
    ) +
    theme_minimal()
  
  # Save plots
  plot_file <- file.path(out_path, paste0("spearman_correlation_N_", num_indiv))
  save_plot(p, plot_file, w = 6, h = 4)
}

# Combine and write results
results_df <- bind_rows(results_list)
print(results_df)

write.csv(results_df, file = file.path(out_path, "spearman_correlation_results.csv"), row.names = FALSE)

# Reproducibility info
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()