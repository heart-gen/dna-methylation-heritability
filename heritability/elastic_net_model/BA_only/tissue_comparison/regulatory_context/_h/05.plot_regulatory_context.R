#### Manuscript-style regulatory-context figures ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
cohort <- ifelse(length(args) >= 1, args[[1]], "BA_only")
population <- ifelse(length(args) >= 2, toupper(args[[2]]), "AA")

base_dir <- here("heritability", "elastic_net_model", cohort,
                 "tissue_comparison", "regulatory_context", "_m")
fig_dir <- file.path(base_dir, "figures", population)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pal_group <- c(
  "Heritable" = "#3B6EA8",
  "Non-heritable" = "#C65146",
  "Low prediction" = "#7A7A7A"
)

theme_regctx <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 1, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size),
      legend.key.height = unit(0.35, "cm"),
      panel.spacing = unit(0.65, "lines")
    )
}

read_many <- function(pattern) {
  files <- list.files(base_dir, pattern = pattern, recursive = TRUE,
                      full.names = TRUE)
  if (length(files) == 0) return(tibble())
  rbindlist(lapply(files, fread), fill = TRUE) |> as_tibble()
}

## Older summaries lacked run_tag (expression implied ABC; PSI implied window_*).
ensure_run_tag <- function(df, default_psi_window_kb = 250L) {
  if (!nrow(df)) return(df)
  if (!"run_tag" %in% names(df)) {
    df$run_tag <- dplyr::case_when(
      df$modality == "expression" ~ "abc",
      df$modality == "psi" ~ paste0("window_", default_psi_window_kb, "kb"),
      TRUE ~ NA_character_
    )
  }
  df
}

regctx_plot_layer <- function(df) {
  df |> mutate(
    regctx_layer = case_when(
      modality == "expression" & run_tag == "abc" ~ "Expression · ABC",
      modality == "expression" &
        grepl("^nearest_gene_window", as.character(run_tag)) ~
        "Expression · nearest gene",
      modality == "psi" ~ "PSI",
      TRUE ~ paste(as.character(modality), as.character(run_tag), sep = " · ")
    )
  )
}

layer_levels <- c(
  "Expression · ABC",
  "Expression · nearest gene",
  "PSI"
)

factor_layer <- function(x) {
  ux <- unique(stats::na.omit(x))
  ord <- c(intersect(layer_levels, ux), setdiff(ux, layer_levels))
  factor(x, levels = ord)
}

save_panel <- function(plot, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = 450)
}

group_df <- read_many("h2_group_association_summary.tsv") |>
  filter(population == !!population) |>
  ensure_run_tag() |>
  regctx_plot_layer() |>
  mutate(regctx_layer = factor_layer(regctx_layer))

bin_df <- read_many("h2_bin_association_summary.tsv") |>
  filter(population == !!population, bin_type == "h2_quintile") |>
  ensure_run_tag() |>
  regctx_plot_layer() |>
  mutate(regctx_layer = factor_layer(regctx_layer))

env_df <- read_many("environment_convergence_enrichment.tsv") |>
  filter(population == !!population) |>
  ensure_run_tag() |>
  regctx_plot_layer() |>
  mutate(regctx_layer = factor_layer(regctx_layer))

prox_group <- read_many("proximity_h2_group_summary.tsv") |>
  filter(population == !!population)
prox_bin <- read_many("proximity_h2_bin_summary.tsv") |>
  filter(population == !!population, bin_type == "h2_quintile")

plots <- list()

if (nrow(group_df) > 0) {
  p_group <- group_df |>
    mutate(
      h2_category = factor(h2_category,
                           levels = c("Heritable", "Non-heritable", "Low prediction"))
    ) |>
    ggplot(aes(h2_category, prop_any_fdr10, fill = h2_category)) +
    geom_col(width = 0.68, color = "white", linewidth = 0.2) +
    facet_grid(regctx_layer ~ tissue) +
    scale_fill_manual(values = pal_group, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "VMRs with linked feature") +
    theme_regctx() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "none")
  plots$group <- p_group
  save_panel(p_group, "panel_h2_group_burden", 6.8, 4.2)
}

if (nrow(bin_df) > 0) {
  pal_layer_line <- c(
    "Expression · ABC" = "#3B6EA8",
    "Expression · nearest gene" = "#7EB8DA",
    "PSI" = "#5A8F61"
  )
  p_bin <- bin_df |>
    ggplot(aes(h2_bin, prop_any_fdr10, group = regctx_layer,
               color = regctx_layer)) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.8) +
    facet_wrap(~ tissue, nrow = 1) +
    scale_color_manual(values = pal_layer_line, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = "h2 quintile", y = "VMRs with linked feature", color = NULL) +
    theme_regctx() +
    theme(legend.position = "top")
  plots$bin <- p_bin
  save_panel(p_bin, "panel_h2_quintile_burden", 6.8, 2.8)
}

if (nrow(env_df) > 0) {
  p_env <- env_df |>
    mutate(
      env_label = gsub("_", " ", env),
      log_or = log2(odds_ratio),
      sig = fdr < 0.10
    ) |>
    ggplot(aes(env_label, tissue, fill = log_or, alpha = sig)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_wrap(~ regctx_layer, nrow = 1) +
    scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#C65146",
                         midpoint = 0, na.value = "grey90",
                         name = "log2 OR") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.45), guide = "none") +
    labs(x = NULL, y = NULL) +
    theme_regctx() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "right")
  plots$env <- p_env
  save_panel(p_env, "panel_environment_convergence", 6.8, 3.6)
}

if (nrow(prox_group) > 0) {
  p_prox_group <- prox_group |>
    mutate(
      h2_category = factor(h2_category,
                           levels = c("Heritable", "Non-heritable", "Low prediction")),
      feature_type = factor(feature_type, levels = c("gene", "psi"))
    ) |>
    ggplot(aes(h2_category, median_distance + 1, color = h2_category,
               group = h2_category)) +
    geom_pointrange(aes(ymin = q25_distance + 1, ymax = q75_distance + 1),
                    position = position_dodge(width = 0.45), linewidth = 0.45) +
    facet_grid(feature_type ~ tissue) +
    scale_color_manual(values = pal_group, drop = FALSE) +
    scale_y_log10(labels = label_number(scale_cut = cut_short_scale())) +
    labs(x = NULL, y = "Nearest distance, bp") +
    theme_regctx() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "none")
  plots$prox_group <- p_prox_group
  save_panel(p_prox_group, "panel_nearest_feature_h2_group", 6.8, 3.4)
}

if (nrow(prox_bin) > 0) {
  p_prox_bin <- prox_bin |>
    mutate(feature_type = factor(feature_type, levels = c("gene", "psi"))) |>
    ggplot(aes(h2_bin, median_distance + 1, group = feature_type,
               color = feature_type)) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.7) +
    facet_wrap(~ tissue, nrow = 1) +
    scale_color_manual(values = c(gene = "#3B6EA8", psi = "#5A8F61")) +
    scale_y_log10(labels = label_number(scale_cut = cut_short_scale())) +
    labs(x = "h2 quintile", y = "Median nearest distance, bp") +
    theme_regctx() +
    theme(legend.position = "top")
  plots$prox_bin <- p_prox_bin
  save_panel(p_prox_bin, "panel_nearest_feature_h2_quintile", 6.8, 2.4)
}

if (length(plots) >= 2) {
  ordered <- plots[c("group", "bin", "prox_group", "prox_bin", "env")]
  ordered <- ordered[!vapply(ordered, is.null, logical(1))]
  main <- wrap_plots(ordered, ncol = 1, guides = "collect") &
    theme(legend.position = "bottom")
  save_panel(main, "main_regulatory_context", 7.2,
             max(4, 2.2 * length(ordered)))
}

message("Saved regulatory-context figures to: ", fig_dir)

#### Reproducibility ####
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
