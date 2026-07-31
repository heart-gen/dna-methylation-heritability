#### Manuscript-ready elastic-net h2 distribution figure across LD decay ####

suppressPackageStartupMessages({
  library(ggplot2)
  library(here)
  library(dplyr)
  library(readr)
})

## Configuration ##
ld_levels <- c("0.8", "0.7", "0.6", "0.5")
category_levels <- c("All sites", "Heritable", "Non-heritable", "Low prediction")
ld_colors <- c(
  "0.8" = "#2B6F8A",
  "0.7" = "#4D8F8A",
  "0.6" = "#B98B3A",
  "0.5" = "#C8624D"
)

## Styling ##
publication_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(
      text = element_text(color = "#1F1F1F"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#E7E7E7", linewidth = 0.3),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "#F5F5F5", color = NA),
      panel.spacing = grid::unit(0.9, "lines"),
      legend.position = "none",
      plot.margin = margin(8, 10, 8, 10)
    )
}

save_plot <- function(p, fn, w, h, dpi = 600) {
  for (ext in c(".png", ".pdf")) {
    ggsave(
      file = paste0(fn, ext), plot = p, width = w, height = h, dpi = dpi
    )
  }
}

## Data helpers ##
get_enet_path <- function(ld_decay) {
  if (ld_decay == "0.8") {
    here(
      "simulation-analysis/elastic-net/sim_200_indiv/_m",
      "simulation_200_summary_elastic-net.tsv"
    )
  } else {
    here(
      "simulation-analysis/elastic-net/ld_decay/_m",
      paste0("simulation_200_", ld_decay, "_summary_elastic-net.tsv")
    )
  }
}

assert_files_exist <- function(paths) {
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths) > 0) {
    stop(
      "Missing required input files:\n",
      paste(missing_paths, collapse = "\n")
    )
  }
}

read_ld_summary <- function(ld_decay) {
  enet_path <- get_enet_path(ld_decay)
  assert_files_exist(enet_path)

  enet <- read_tsv(enet_path, show_col_types = FALSE)

  if (!"LD_Decay" %in% names(enet)) {
    enet <- enet |> mutate(LD_Decay = ld_decay)
  }

  enet |>
    mutate(LD_Decay = factor(as.character(LD_Decay), levels = ld_levels, ordered = TRUE))
}

categorize_sites <- function(enet) {
  bind_rows(
    enet |> mutate(h2_category = "All sites"),
    enet |>
      filter(h2_unscaled >= 0.1, r_squared_cv > 0.75) |>
      mutate(h2_category = "Heritable"),
    enet |>
      filter(h2_unscaled < 0.1, r_squared_cv > 0.75) |>
      mutate(h2_category = "Non-heritable"),
    enet |>
      filter(!is.finite(r_squared_cv) | r_squared_cv <= 0.75) |>
      mutate(h2_category = "Low prediction")
  ) |>
    mutate(h2_category = factor(h2_category, levels = category_levels))
}

## Main ##
out_path <- here("simulation-analysis/comparison/h2_distribution/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

assert_files_exist(vapply(ld_levels, get_enet_path, FUN.VALUE = character(1)))

plot_data <- bind_rows(lapply(ld_levels, function(ld_decay) {
  categorize_sites(read_ld_summary(ld_decay))
})) |>
  mutate(
    LD_Decay = factor(LD_Decay, levels = ld_levels, ordered = TRUE),
    h2_category = factor(h2_category, levels = category_levels)
  ) |>
  filter(is.finite(h2_unscaled))

summary_df <- plot_data |>
  group_by(LD_Decay, h2_category) |>
  summarise(
    n = n(), mean_h2 = mean(h2_unscaled, na.rm = TRUE),
    median_h2 = median(h2_unscaled, na.rm = TRUE), .groups = "drop"
  ) |> arrange(LD_Decay, h2_category)

y_upper <- max(plot_data$h2_unscaled, na.rm = TRUE) * 1.03

p <- ggplot(plot_data, aes(x = LD_Decay, y = h2_unscaled, fill = LD_Decay)) +
  geom_hline(
    yintercept = 0.1, color = "#7A7A7A", linetype = "22", linewidth = 0.35
  ) +
  geom_violin(
    trim = FALSE, scale = "width", linewidth = 0.25, color = "#F7F7F7", alpha = 0.95
  ) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", 
               color = "#2A2A2A", linewidth = 0.32, show.legend = FALSE
  ) +
  stat_summary(
    fun = median, geom = "point", size = 1.8, stroke = 0.3, 
    fill = "white", color = "#2A2A2A"
  ) +
  facet_wrap(~h2_category, ncol = 2) +
  scale_fill_manual(values = ld_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  coord_cartesian(ylim = c(0, y_upper)) +
  labs(x = "LD decay", y = expression("Estimated " * h^2)) +
  publication_theme()

fn_hist <- file.path(out_path, "enet_h2_distribution_ld_decay_200")
save_plot(p, fn_hist, w = 9.4, h = 7.2)

print(summary_df)

write_csv(summary_df, file.path(out_path, "h2_summary_ld_decay_200.csv"))

## Reproducibility ##
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
