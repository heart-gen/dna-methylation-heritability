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
out_path <- here("simulation-analysis/comparison/summary_counts/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# Read in and filter summary tables
num_indivs <- c(100, 150, 200, 250, 500, 1000, 5000)
enet_list <- list()

for (num_indiv in num_indivs) {
  enet_file <- here("simulation-analysis/elastic-net", paste0("sim_", num_indiv, "_indiv"), "_m",
                    paste0("simulation_", num_indiv, "_summary_elastic-net.tsv"))
  enet <- read.table(enet_file, sep = "\t", header = TRUE)
  enet <- na.omit(enet)
  enet$N <- num_indiv
  enet_list[[as.character(num_indiv)]] <- enet
}

enet_all <- bind_rows(enet_list)

enet_all <- enet_all %>%
  mutate(method = "Elastic Net", N = as.character(N),
         h2_category = case_when(
           r_squared_cv <= 0.75 ~ "Low prediction",
           h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
           h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
         ))

gcta_path <- here("simulation-analysis/gcta", "_m", "summary", "greml_summary.tsv")
gcta <- read.table(gcta_path, sep = "\t", header = TRUE)
gcta <- gcta %>%
  mutate(method = "GREML-LDMS", N = as.character(N),
         p_adjusted_fdr = p.adjust(Pval_Variance, method = "fdr"),
         h2_category = case_when(
           p_adjusted_fdr >= 0.05 ~ "Low prediction",
           Sum.of.V.G._Vp_Variance < 0.1 & p_adjusted_fdr < 0.05 ~ "Non-heritable",
           Sum.of.V.G._Vp_Variance > 0.1 & p_adjusted_fdr < 0.05 ~ "Heritable"
         ))


target_path <- here("inputs", "simulated-data", "_m", paste0("sim_", num_indiv, "_indiv"), 
                    "snp_phenotype_mapping.tsv")
target <- read.table(target_path, sep = "\t", header = TRUE)
target <- target %>%
  mutate(N = "Target", method = "Target",
         h2_category = case_when(
           target_heritability < 0.1 ~ "Non-heritable",
           target_heritability >= 0.1 ~ "Heritable"
         ))

# Summarize categories for each metho
gcta_summary <- gcta %>%
  group_by(N, method, h2_category) %>%
  summarise(count = n(), .groups = "drop")

enet_summary <- enet_all %>%
  group_by(N, method, h2_category) %>%
  summarise(count = n(), .groups = "drop")

target_summary <- target %>%
  group_by(N, method, h2_category) %>%
  summarise(count = n(), .groups = "drop")

# Combine summary counts
all_summary <- bind_rows(gcta_summary, enet_summary, target_summary)

# Get failed instances (for gcta)
non_failed_summary <- all_summary %>%
  filter(h2_category != "Failed") %>%
  group_by(N, method) %>%
  summarise(non_failed_total = sum(count), .groups = "drop")

failed_summary <- non_failed_summary %>%
  mutate(h2_category = "Failed", count = 1000 - non_failed_total) %>%
  select(N, method, h2_category, count)

gcta_n100_failed <- tibble(
  N = "100", method = "GREML-LDMS",
  h2_category = "Failed", count = 1000
)

# Add failed instances to summary
all_summary_with_failed <- bind_rows(all_summary, failed_summary, gcta_n100_failed)

# Order correctly for plotting
all_summary_with_failed <- all_summary_with_failed %>%
  mutate(
    h2_category = factor(h2_category, 
                         levels = c("Heritable", "Non-heritable", "Low prediction", "Failed")),
    method = factor(method, levels = c("Target", "Elastic Net", "GREML-LDMS")),
    N = factor(N, levels = c("Target", as.character(sort(num_indivs))))
  ) %>%
  filter(!is.na(N)) 

# Plot grouped stacked bar plot
p <- ggplot(all_summary_with_failed, aes(x = N, y = count, fill = h2_category)) +
  geom_bar(
    aes(group = method), 
    stat = "identity",
    position = position_stack(), 
    color = "black"
  ) +
  facet_grid(. ~ method, switch = "x", scales = "free_x", space = "free_x") +
  scale_fill_manual(
    values = c("Heritable" = "#497C8A",
               "Non-heritable" = "#8CA77B",
               "Low prediction" = "#E3A27F",
               "Failed" = "#B34040")
  ) +
  labs(
    x = "Sample Size (N)",
    y = "Count",
    fill = "Category"
  ) +
  theme_minimal(base_size = 20) +
  #font("xy.title", face = "bold", size = 14) + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.title = element_text(hjust = 0.5),
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save plot
plot_file <- file.path(out_path, "simulated_data_stacked")
save_plot(p, plot_file, w = 10, h = 6, dpi = 300)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()