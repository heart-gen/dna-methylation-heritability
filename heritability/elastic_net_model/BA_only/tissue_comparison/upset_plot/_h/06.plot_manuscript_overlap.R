## This script assembles the manuscript-ready F0.25 overlap figure

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

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
      legend.position = "none"
    )
}

region_from_membership <- function(region_key) {
  dplyr::case_when(
    region_key == "includes_caudate" ~ "Caudate",
    region_key == "includes_dlpfc" ~ "DLPFC",
    region_key == "includes_hippocampus" ~ "Hippocampus",
    TRUE ~ region_key
  )
}

build_panel_a <- function(set_counts) {
  h2_levels <- c("All VMRs", "Heritable", "Non-heritable", "Low prediction")
  panel_metrics <- set_counts %>%
    group_by(h2_label) %>%
    summarise(max_count = max(count), .groups = "drop") %>%
    mutate(
      h2_label = factor(h2_label, levels = h2_levels),
      y_caudate = -0.08 * max_count,
      y_dlpfc = -0.17 * max_count,
      y_hippocampus = -0.26 * max_count,
      y_bottom = -0.35 * max_count
    )

  bar_df <- set_counts %>%
    left_join(panel_metrics, by = "h2_label") %>%
    mutate(
      h2_label = factor(h2_label, levels = h2_levels),
      label_y = count + 0.03 * max_count
    )

  matrix_df <- set_counts %>%
    left_join(panel_metrics, by = "h2_label") %>%
    select(
      h2_label,
      set_order,
      includes_caudate,
      includes_dlpfc,
      includes_hippocampus,
      y_caudate,
      y_dlpfc,
      y_hippocampus
    ) %>%
    pivot_longer(
      cols = starts_with("includes_"),
      names_to = "region_key",
      values_to = "included"
    ) %>%
    mutate(
      h2_label = factor(h2_label, levels = h2_levels),
      region_label = region_from_membership(region_key),
      region_label = factor(region_label, levels = c("Caudate", "DLPFC", "Hippocampus")),
      y = case_when(
        region_label == "Caudate" ~ y_caudate,
        region_label == "DLPFC" ~ y_dlpfc,
        TRUE ~ y_hippocampus
      ),
      included_flag = ifelse(included == 1, "Included", "Not included")
    )

  line_df <- matrix_df %>%
    filter(included == 1) %>%
    group_by(h2_label, set_order) %>%
    summarise(
      ymin = min(y),
      ymax = max(y),
      .groups = "drop"
    ) %>%
    filter(ymax > ymin)

  all_panel_metrics <- panel_metrics %>%
    filter(h2_label == "All VMRs")

  label_df <- data.frame(
    h2_label = factor(rep("All VMRs", 3), levels = h2_levels),
    set_order = rep(0.1, 3),
    label = c("Caudate", "DLPFC", "Hippocampus"),
    y = c(all_panel_metrics$y_caudate, all_panel_metrics$y_dlpfc, all_panel_metrics$y_hippocampus)
  )

  ggplot(bar_df, aes(x = set_order, y = count)) +
    geom_col(width = 0.72, fill = "#2F5964") +
    geom_text(aes(y = label_y, label = comma(count)), size = 3) +
    geom_segment(
      data = line_df,
      aes(x = set_order, xend = set_order, y = ymin, yend = ymax),
      inherit.aes = FALSE,
      color = "#76888A",
      linewidth = 0.5
    ) +
    geom_point(
      data = matrix_df,
      aes(x = set_order, y = y, fill = included_flag),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.8,
      stroke = 0.4,
      color = "#5E6667"
    ) +
    geom_text(
      data = label_df,
      aes(x = set_order, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      size = 4
    ) +
    facet_wrap(~h2_label, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c("Included" = "#2F5964", "Not included" = "#F4F1EA")) +
    scale_x_continuous(
      breaks = 1:7,
      labels = rep("", 7),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_y_continuous(
      labels = function(x) ifelse(x < 0, "", label_number(big.mark = ",")(x)),
      expand = expansion(mult = c(0.08, 0.16))
    ) +
    labs(x = NULL, y = "Shared VMR count") +
    theme_overlap(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(5, "mm"),
      plot.margin = margin(5.5, 12, 5.5, 38)
    ) +
    coord_cartesian(clip = "off")
}

build_heatmap <- function(df, value_col, label_col, fill_scale) {
  h2_levels <- c("All VMRs", "Heritable", "Non-heritable", "Low prediction")
  pair_levels <- c("Caudate\nDLPFC", "Caudate\nHippocampus", "Hippocampus\nDLPFC")

  plot_df <- df %>%
    mutate(
      h2_label = factor(h2_label, levels = rev(h2_levels)),
      pair_label_multiline = factor(pair_label_multiline, levels = pair_levels)
    )

  ggplot(plot_df, aes(x = pair_label_multiline, y = h2_label, fill = .data[[value_col]])) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = .data[[label_col]]), size = 3.1, lineheight = 0.96) +
    fill_scale +
    labs(x = NULL, y = NULL) +
    theme_overlap(base_size = 11) +
    theme(
      axis.text.x = element_text(face = "bold", lineheight = 0.95),
      axis.text.y = element_text(face = "bold"),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
}

## Main
out_dir <- here("heritability/elastic_net_model/BA_only/tissue_comparison/upset_plot/_m/manuscript_overlap_figure")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

set_counts <- fread(file.path(out_dir, "f0.25_set_counts.tsv"))
reciprocal_summary <- fread(file.path(out_dir, "f0.25_reciprocal_overlap_summary.tsv"))
h2_summary <- fread(file.path(out_dir, "f0.25_h2_concordance_summary.tsv"))

reciprocal_summary <- reciprocal_summary %>%
  mutate(
    heatmap_label = paste0(number(median_reciprocal, accuracy = 0.1), "%\n", "n=", comma(n_pairs))
  )

h2_summary <- h2_summary %>%
  mutate(
    heatmap_label = paste0(sprintf("%.2f", spearman_rho), "\n", "n=", comma(n_pairs))
  )

panel_a <- build_panel_a(set_counts)

panel_b <- build_heatmap(
  reciprocal_summary,
  value_col = "median_reciprocal",
  label_col = "heatmap_label",
  fill_scale = scale_fill_gradientn(
    colours = c("#F4F1EA", "#D7E3DC", "#8FAFAC", "#2F5964"),
    limits = c(25, 100),
    oob = squish
  )
)

panel_c <- build_heatmap(
  h2_summary,
  value_col = "spearman_rho",
  label_col = "heatmap_label",
  fill_scale = scale_fill_gradient2(
    low = "#D6E0E3",
    mid = "#F7F4EF",
    high = "#B8623F",
    midpoint = 0,
    limits = c(-0.05, 0.8),
    oob = squish
  )
)

assembled_figure <- wrap_plots(
  A = wrap_elements(full = panel_a),
  B = panel_b,
  C = panel_c,
  design = "
AA
BC
",
  heights = c(1.55, 1)
) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14, color = "#1A1A1A")
    )
  )

save_plot(panel_a, file.path(out_dir, "panel_A_f0.25_upset_summary"), 12, 5.3)
save_plot(panel_b, file.path(out_dir, "panel_B_f0.25_reciprocal_overlap"), 5.5, 4.6)
save_plot(panel_c, file.path(out_dir, "panel_C_f0.25_h2_concordance"), 5.5, 4.6)
save_plot(assembled_figure, file.path(out_dir, "f0.25_overlap_manuscript_figure"), 12, 8.6)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
if (requireNamespace("sessioninfo", quietly = TRUE)) {
  sessioninfo::session_info()
}
