## This script plots concordance for quintiles of matched VMRs in AA and EA

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## Function
save_plot <- function(plot_obj, file_base, width, height, dpi = 400) {
  for (ext in c(".pdf", ".png")) {
    ggsave(
      filename = paste0(file_base, ext),
      plot = plot_obj,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }
}

theme_overlap <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "#1A1A1A"),
      strip.text = element_text(face = "bold", color = "#1A1A1A"),
      strip.background = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

build_heatmap <- function(df, fill_scale) {

  ggplot(df, aes(x = h2_quintile_AA, y = h2_quintile_EA, fill = n)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = comma(n)), size = 3.1, lineheight = 0.96) +
    fill_scale +
    facet_wrap(~ tissue, nrow = 1) +
    labs(x = "BA h2 quintile", y = "WA h2 quintile") +
    theme_overlap(base_size = 11) +
    theme(
      axis.text.x = element_text(face = "bold", lineheight = 0.95),
      axis.text.y = element_text(face = "bold"),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

## Main
out_path <- here("heritability/elastic_net_model/all_individuals/tissue_comparison/quintile_summary/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

tissues <- c("caudate", "dlpfc", "hippocampus")

quintile_all <- list()

for (tissue in tissues){

    quintile_summary <- fread(file.path(out_path, paste0("quintile_match_summary_AA_EA_", tissue, ".csv"))) %>%
    mutate(tissue = tissue)


    quintile_all[[tissue]] <- quintile_summary
}

quintile_summary_all <- bind_rows(quintile_all)

quintile_levels <- c("Q1", "Q2", "Q3", "Q4", "Q5")

quintile_summary_all <- quintile_summary_all %>%
  mutate(
    h2_quintile_AA = factor(h2_quintile_AA, levels = quintile_levels),
    h2_quintile_EA = factor(h2_quintile_EA, levels = quintile_levels),
    tissue = factor(
      tissue,
      levels = c("caudate", "dlpfc", "hippocampus"),
      labels = c("Caudate", "DLPFC", "Hippocampus")
    )
  )

heatmap <- build_heatmap(
  quintile_summary_all,
  fill_scale = scale_fill_gradientn(
    colours = c("#F4F1EA", "#D7E3DC", "#8FAFAC", "#2F5964"),
    limits = c(0, max(quintile_summary_all$n)),
    oob = squish
  )
)

save_plot(heatmap, file.path(out_path, "quintile_match_summary"), 6, 4)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()