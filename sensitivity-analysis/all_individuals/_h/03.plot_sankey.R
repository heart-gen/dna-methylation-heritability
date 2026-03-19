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
out_path <- here("sensitivity-analysis/all_individuals/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# Read in and filter summary tables
tissues <- c("caudate", "hippocampus", "dlpfc")
populations <- c("AA", "EA")

enet_list <- list()

for (tissue in tissues) {
  for (pop in populations) {
    enet_file <- here("heritability/elastic_net_model/all_individuals", tissue, "_m",
                      paste0(tissue, "_summary_elastic-net_", pop, ".tsv"))
    enet <- read.table(enet_file, sep = "\t", header = TRUE)
    enet <- na.omit(enet)
    enet_list[[paste(tissue, pop, sep = "_")]] <- enet
  }
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

nh_aa_summary <- all_summary %>%
  filter(race == "AA", h2_category == "Non-heritable") %>%
  select(chrom, start, end, region, r2_threshold, h2_category, race)

aa_to_ea <- nh_aa_summary %>%
  left_join(all_summary %>% filter(race == "EA") %>% 
              select(chrom, start, end, region, r2_threshold, h2_category, race),
            by = c("region", "r2_threshold", "chrom", "start", "end"),
            suffix = c("_AA", "_EA")) %>%
  filter(!is.na(h2_category_EA))

summary <- aa_to_ea %>%
  group_by(region, r2_threshold, h2_category_AA, h2_category_EA) %>%
  summarise(count = n(), .groups = "drop")

p_sankey <- ggplot(summary, aes(axis1 = h2_category_AA, 
                                axis2 = factor(region),
                                axis3 = h2_category_EA, y = count)) +
  geom_alluvium(aes(fill = h2_category_EA), width = 1/12) + 
  geom_stratum(width = 1/8, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = paste0(after_stat(stratum), "\n", after_stat(count))), 
            angle = 45, size = 4, nudge_y = 0.05) +
  facet_grid(~ r2_threshold, scales = "free_y") +
  scale_fill_manual(
    values = c("Heritable" = "#497C8A",
               "Non-heritable" = "#8CA77B",
               "Low prediction" = "#E3A27F")
  ) +
  labs(
    x = NULL,
    y = "Count",
    fill = "EA Category",
    title = "AA Non-heritable → EA Categories"
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

print(p_sankey)

# Save plot
plot_file <- file.path(out_path, "AA_to_EA_sankey")
save_plot(p_sankey, plot_file, w = 12, h = 6, dpi = 300)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()