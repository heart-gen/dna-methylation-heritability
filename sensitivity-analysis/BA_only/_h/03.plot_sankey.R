# Track heritability category classification across AA and EA

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggalluvial)
})

## Function
save_plot <- function(p, fn, w, h, dpi){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn,ext), plot=p, width=w, height=h, dpi=dpi)
  }
}

## Main 
# Create output directory
out_path <- here("sensitivity-analysis/BA_only/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# Read in and filter summary tables
tissues <- c("caudate", "hippocampus", "dlpfc")

enet_list <- list()

for (tissue in tissues) {
    enet_file <- here("heritability/elastic_net_model/all_individuals", tissue, "_m",
                      paste0(tissue, "_summary_elastic-net.tsv"))
    enet <- read.table(enet_file, sep = "\t", header = TRUE)
    enet <- na.omit(enet)
    enet_list[[paste(tissue, pop, sep = "_")]] <- enet
}

enet_all <- bind_rows(enet_list)

# Get summary for each r2 threshold
enet_high <- enet_all %>%
  mutate(h2_category = case_when(
    r_squared_cv <= 0.75 ~ "Low prediction",
    h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
    h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
  ),
  r2_threshold = 0.75)

enet_high <- enet_high %>%
  mutate(h2_category = factor(h2_category,
                              levels = c("Heritable", "Non-heritable", "Low prediction")),
         region = recode(region, "caudate" = "Caudate"))

enet_low <- enet_all %>%
  mutate(h2_category = case_when(
    r_squared_cv <= 0.3 ~ "Low prediction",
    h2_unscaled < 0.1 & r_squared_cv > 0.3 ~ "Non-heritable",
    h2_unscaled >= 0.1 & r_squared_cv > 0.3 ~ "Heritable"
  ),
  r2_threshold = 0.3)

enet_low <- enet_low %>%
  mutate(h2_category = factor(h2_category,
                              levels = c("Heritable", "Non-heritable", "Low prediction")),
         region = recode(region, "caudate" = "Caudate"))

# Summarize by r2 threshold
all_summary <- bind_rows(enet_high, enet_low)

low_pred_summary <- all_summary %>%
  pivot_wider(
    names_from = r2_threshold,
    values_from = h2_category,
    names_prefix = "r2_"
  ) %>%
  filter(r2_0.75 == "Low prediction")

summary <- low_pred_summary %>%
  group_by(region, r2_0.75, r2_0.3) %>%
  summarise(count = n(), .groups = "drop")

p_sankey <- ggplot(summary, aes(axis1 = r2_0.75, 
                                axis2 = factor(region),
                                axis3 = r2_0.3, y = count)) +
  geom_alluvium(aes(fill = r2_0.3), width = 1/12) + 
  geom_stratum(width = 1/3, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", 
            aes(label = paste0(after_stat(stratum), "\n", 
                               after_stat(count), "\n", 
                               scales::percent(after_stat(count / sum(count)), accuracy = 0.01))), 
            size = 4, nudge_y = 0.05) +
  scale_fill_manual(
    values = c("Heritable" = "#497C8A",
               "Non-heritable" = "#8CA77B",
               "Low prediction" = "#E3A27F")
  ) +
  labs(
    x = NULL,
    y = "Count",
    fill = "0.3 Category",
    title = "0.75 Low Prediction → 0.3 Categories"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.title = element_text(hjust = 0.5),
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
  )

# Save plot
plot_file <- file.path(out_path, "0.75_to_0.3_sankey")
save_plot(p_sankey, plot_file, w = 12, h = 6, dpi = 300)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()