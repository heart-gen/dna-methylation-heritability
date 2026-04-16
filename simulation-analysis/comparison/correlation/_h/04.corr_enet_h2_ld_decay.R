##### Elastic-net LD-decay correlation figure: boosting_hybrid vs. joint_ridge #####
suppressPackageStartupMessages({
  library(ggplot2)
  library(here)
  library(dplyr)
  library(readr)
})

## Configuration ##
ld_levels       <- c("0.8", "0.7", "0.6", "0.5")
enet_methods    <- c("boosting_hybrid", "joint_ridge")
method_labels   <- c(boosting_hybrid = "Boosting-hybrid", joint_ridge = "Joint-ridge")
category_levels <- c("Heritable", "Non-heritable", "Low prediction")
category_colors <- c(
  "Heritable"      = "#2B6F8A",
  "Non-heritable"  = "#7E9B63",
  "Low prediction" = "#D38A5C"
)

## Styling ##
publication_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(
      text = element_text(color = "#1F1F1F"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E7E7E7", linewidth = 0.3),
      strip.text = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "#F5F5F5", color = NA),
      strip.placement = "outside",
      panel.spacing = grid::unit(0.75, "lines"),
      legend.position = "none",
      plot.margin = margin(8, 10, 8, 10)
    )
}

save_plot <- function(p, fn, w, h, dpi = 600) {
  for (ext in c(".png", ".pdf")) {
    ggsave(
      file = paste0(fn, ext), plot = p, width = w,
      height = h, dpi = dpi
    )
  }
}

## Data helpers ##
get_enet_path <- function(ld_decay, enet_method) {
  if (ld_decay == "0.8") {
    here(
      "simulation-analysis/elastic-net/sim_200_indiv/_m",
      "simulation_200_summary_elastic-net.tsv"
    )
  } else {
    here(
      "simulation-analysis/elastic-net/ld_decay/_m",
      paste0("simulation_200_", ld_decay, "_", enet_method,
             "_summary_elastic-net.tsv")
    )
  }
}

get_target_path <- function(ld_decay) {
  if (ld_decay == "0.8") {
    here(
      "inputs/simulated-data/_m/sim_200_indiv",
      "snp_phenotype_mapping.tsv"
    )
  } else {
    here(
      "inputs/simulated-data/_m",
      paste0("ld_", ld_decay, "_sim_200_indiv"),
      "snp_phenotype_mapping.tsv"
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

read_ld_summary <- function(ld_decay, enet_method) {
  enet_path   <- get_enet_path(ld_decay, enet_method)
  target_path <- get_target_path(ld_decay)
  assert_files_exist(c(enet_path, target_path))

  enet   <- read_tsv(enet_path,   show_col_types = FALSE)
  target <- read_tsv(target_path, show_col_types = FALSE)

  if (!"LD_Decay" %in% names(enet)) {
    enet <- enet |> mutate(LD_Decay = ld_decay)
  }

  merged <- inner_join(enet, target, by = c("pheno_id" = "phenotype_id"))

  if (nrow(merged) != nrow(enet)) {
    stop(
      "Elastic-net and target tables did not join cleanly for LD decay ",
      ld_decay, " method ", enet_method,
      ". Expected ", nrow(enet), " rows but found ", nrow(merged), "."
    )
  }

  merged |>
    mutate(
      LD_Decay    = factor(as.character(LD_Decay), levels = ld_levels, ordered = TRUE),
      enet_method = enet_method,
      method_label = method_labels[enet_method],
      h2_category = case_when(
        !is.finite(r_squared_cv)                 ~ "Low prediction",
        r_squared_cv <= 0.75                     ~ "Low prediction",
        target_heritability < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
        target_heritability >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
      ),
      h2_category = factor(h2_category, levels = category_levels)
    )
}

get_spearman_stats <- function(df) {
  usable_df <- df |>
    filter(is.finite(target_heritability), is.finite(h2_unscaled))

  if (
    nrow(usable_df) < 3 ||
      n_distinct(usable_df$target_heritability) < 2 ||
      n_distinct(usable_df$h2_unscaled) < 2
  ) {
    return(tibble(spearman_rho = NA_real_, p_value = NA_real_, n = nrow(usable_df)))
  }

  spearman_test <- cor.test(
    usable_df$target_heritability, usable_df$h2_unscaled,
    method = "spearman"
  )
  tibble(
    spearman_rho = unname(spearman_test$estimate),
    p_value      = spearman_test$p.value,
    n            = nrow(usable_df)
  )
}

format_stat_label <- function(rho, n) {
  if (is.na(rho)) {
    paste0("n = ", n)
  } else {
    paste0("rho = ", formatC(rho, format = "f", digits = 2), "\n", "n = ", n)
  }
}

format_ld_strip <- function(x) paste0("LD decay = ", x)

## Main ##
out_path <- here("simulation-analysis/comparison/correlation/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

all_required_paths <- unlist(lapply(enet_methods, function(m) {
  unlist(lapply(ld_levels, function(ld) c(get_enet_path(ld, m), get_target_path(ld))))
}))
assert_files_exist(all_required_paths)

all_data <- bind_rows(lapply(enet_methods, function(m) {
  bind_rows(lapply(ld_levels, function(ld) read_ld_summary(ld, m)))
})) |>
  mutate(
    LD_Decay     = factor(LD_Decay, levels = ld_levels, ordered = TRUE),
    method_label = factor(method_labels[enet_method],
                          levels = unname(method_labels)),
    h2_category  = factor(h2_category, levels = category_levels)
  )

plot_data <- all_data |>
  filter(is.finite(target_heritability), is.finite(h2_unscaled))

results_df <- all_data |>
  group_by(LD_Decay, method_label, h2_category) |>
  group_modify(~get_spearman_stats(.x)) |>
  ungroup() |>
  mutate(
    LD_Decay     = factor(LD_Decay, levels = ld_levels, ordered = TRUE),
    method_label = factor(method_label, levels = unname(method_labels)),
    h2_category  = factor(h2_category, levels = category_levels)
  ) |>
  arrange(LD_Decay, method_label, h2_category)

x_upper <- max(plot_data$target_heritability, na.rm = TRUE) * 1.03
y_upper <- max(plot_data$h2_unscaled,         na.rm = TRUE) * 1.03

annotation_df <- results_df |>
  mutate(
    label   = mapply(format_stat_label, spearman_rho, n),
    label_x = x_upper * 0.04,
    label_y = y_upper * 0.96
  )

p <- ggplot(plot_data, aes(x = target_heritability, y = h2_unscaled)) +
  geom_abline(
    intercept = 0, slope = 1, color = "#C8C8C8", linetype = "dotted",
    linewidth = 0.35
  ) +
  geom_vline(xintercept = 0.1, color = "#7A7A7A", linetype = "22", linewidth = 0.35) +
  geom_hline(yintercept = 0.1, color = "#7A7A7A", linetype = "22", linewidth = 0.35) +
  geom_point(
    aes(color = h2_category), alpha = 0.58, size = 0.85, show.legend = FALSE
  ) +
  geom_smooth(
    aes(color = h2_category), method = "lm", formula = y ~ x,
    se = FALSE, linewidth = 0.6, show.legend = FALSE
  ) +
  geom_text(
    data = annotation_df, aes(x = label_x, y = label_y, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.0,
    lineheight = 0.95, color = "#2A2A2A"
  ) +
  facet_grid(
    rows = vars(LD_Decay), cols = vars(h2_category),
    switch = "y",
    labeller = labeller(
      LD_Decay    = format_ld_strip,
      h2_category = label_value
    )
  ) +
  scale_color_manual(values = category_colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, x_upper), ylim = c(0, y_upper), clip = "off") +
  labs(
    x = expression("True " * h^2),
    y = expression("Estimated " * h^2)
  ) +
  publication_theme()

## One panel per method
for (m in enet_methods) {
  p_m <- p %+% filter(plot_data, enet_method == m) +
    labs(subtitle = method_labels[m])
  ann_m <- filter(annotation_df, method_label == method_labels[m])
  p_m <- p_m + geom_text(
    data = ann_m, aes(x = label_x, y = label_y, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.0,
    lineheight = 0.95, color = "#2A2A2A"
  )
  plot_file_m <- file.path(
    out_path,
    paste0("elastic_net_ld_decay_correlation_", m)
  )
  save_plot(p_m, plot_file_m, w = 10.6, h = 10.2)
}

## Combined panel (facet by method_label × h2_category)
p_combined <- ggplot(plot_data, aes(x = target_heritability, y = h2_unscaled)) +
  geom_abline(
    intercept = 0, slope = 1, color = "#C8C8C8", linetype = "dotted",
    linewidth = 0.35
  ) +
  geom_vline(xintercept = 0.1, color = "#7A7A7A", linetype = "22", linewidth = 0.35) +
  geom_hline(yintercept = 0.1, color = "#7A7A7A", linetype = "22", linewidth = 0.35) +
  geom_point(
    aes(color = h2_category), alpha = 0.45, size = 0.6, show.legend = FALSE
  ) +
  geom_smooth(
    aes(color = h2_category), method = "lm", formula = y ~ x,
    se = FALSE, linewidth = 0.5, show.legend = FALSE
  ) +
  geom_text(
    data = annotation_df, aes(x = label_x, y = label_y, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.5,
    lineheight = 0.9, color = "#2A2A2A"
  ) +
  facet_grid(
    rows = vars(LD_Decay), cols = vars(method_label),
    switch = "y",
    labeller = labeller(LD_Decay = format_ld_strip, method_label = label_value)
  ) +
  scale_color_manual(values = category_colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, x_upper), ylim = c(0, y_upper), clip = "off") +
  labs(
    x = expression("True " * h^2),
    y = expression("Estimated " * h^2)
  ) +
  publication_theme()

plot_file <- file.path(out_path, "elastic_net_ld_decay_correlation_combined")
save_plot(p_combined, plot_file, w = 14.0, h = 10.2)

print(results_df)

write_csv(
  results_df,
  file.path(out_path, "elastic_net_ld_decay_spearman_correlation_results.csv")
)

## Reproducibility ##
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
