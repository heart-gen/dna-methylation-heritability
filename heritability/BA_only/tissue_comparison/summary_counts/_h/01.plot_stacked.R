# Summarize heritability classifications using stacked barplot

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
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
out_path <- here("heritability/elastic_net_model/BA_only/tissue_comparison/summary_counts/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# Read in and filter summary tables
tissues <- c("caudate", "hippocampus", "dlpfc")
enet_list <- list()

for (tissue in tissues) {
  enet_file <- here("heritability/elastic_net_model/BA_only", tissue, "_m",
                    paste0(tissue, "_summary_elastic-net.tsv"))
  enet <- read.table(enet_file, sep = "\t", header = TRUE)
  enet <- na.omit(enet)
  enet_list[[as.character(tissue)]] <- enet
}

enet_all <- bind_rows(enet_list)

enet_all <- enet_all %>%
  mutate(h2_category = case_when(
           r_squared_cv <= 0.3 ~ "Low prediction",
           h2_unscaled < 0.1 & r_squared_cv > 0.3 ~ "Non-heritable",
           h2_unscaled >= 0.1 & r_squared_cv > 0.3 ~ "Heritable"
         ))

enet_all <- enet_all %>%
  mutate(h2_category = factor(h2_category,
                              levels = c("Heritable", "Non-heritable", "Low prediction")),
         region = recode(region, "caudate" = "Caudate"))

enet_summary <- enet_all %>%
  group_by(region, h2_category) %>%
  summarise(count = n(), .groups = "drop")
print(enet_summary)

# Plot grouped stacked bar plot
p <- ggplot(enet_summary, aes(x = region, y = count, fill = h2_category)) +
  geom_bar(
    aes(group = region), 
    stat = "identity",
    position = position_stack(), 
    color = "black"
  ) +
  scale_fill_manual(
    values = c("Heritable" = "#497C8A",
               "Non-heritable" = "#8CA77B",
               "Low prediction" = "#E3A27F")
  ) +
  labs(
    x = "Brain Region",
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
plot_file <- file.path(out_path, "vmr_summary_stacked")
save_plot(p, plot_file, w = 6, h = 5, dpi = 300)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()