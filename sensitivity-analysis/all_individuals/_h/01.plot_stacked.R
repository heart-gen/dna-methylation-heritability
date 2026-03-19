# Summarize heritability classifications using stacked barplot

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
})

## Function
save_plot <- function(p, fn, w, h, dpi){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn,ext), plot=p, width=w, height=h, dpi=dpi)
  }
}

plot_stacked <- function(enet_summary, population){
  enet_summary <- enet_summary %>%
    filter(race == population)
  
  p <- ggplot(enet_summary, aes(x = region, y = count, fill = h2_category)) +
    geom_bar(
      aes(group = interaction(region, race)), 
      stat = "identity",
      position = position_stack(), 
      color = "black"
    ) +
    facet_grid(. ~ r2_threshold, switch = "x", scales = "free_x", space = "free_x") +
    scale_fill_manual(
      values = c("Heritable" = "#497C8A",
                 "Non-heritable" = "#8CA77B",
                 "Low prediction" = "#E3A27F")
    ) +
    labs(
      x = "Brain Region",
      y = "Count",
      fill = "Category",
      title = population
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

# Summarize number of lost VMRs
loss <- enet_all %>%
  mutate(
    low = r_squared_cv >= 0.3,
    high = r_squared_cv >= 0.75,
    lost = low & !high
  )

loss_summary <- loss %>%
  group_by(race, region) %>%
  summarise(count = sum(lost), .groups = "drop")

# Write summary counts to file
write.csv(loss_summary, file = file.path(out_path, "summary_loss_count.csv"), 
          row.names = FALSE)

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
enet_summary <- all_summary %>%
  group_by(race, r2_threshold, region, h2_category) %>%
  summarise(count = n(), .groups = "drop")

# Write summary counts to file
write.csv(enet_summary, file = file.path(out_path, "summary_count.csv"), 
          row.names = FALSE)

# Plot grouped stacked bar plots
p_AA <- plot_stacked(enet_summary, "AA")
p_EA <- plot_stacked(enet_summary, "EA")
p <- p_AA + p_EA +
  plot_layout(guides = "collect") + 
  plot_annotation(theme = theme(legend.position = "right"))

# Save plot
plot_file <- file.path(out_path, "vmr_summary_stacked")
save_plot(p, plot_file, w = 12, h = 6, dpi = 300)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()