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

source(here("heritability", "elastic_net_model", "all_individuals",
            "tissue_comparison", "regulatory_context", "_h",
            "00.regulatory_context_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
cohort <- ifelse(length(args) >= 1, args[[1]], "all_individuals")
population <- ifelse(length(args) >= 2, toupper(args[[2]]), "AA")
vmr_set <- ifelse(length(args) >= 3, args[[3]], "shared")
vmr_set <- validate_vmr_set(cohort, population, vmr_set)
plot_populations <- shared_vmr_plot_populations(population, vmr_set)
canonical_pop <- if (normalize_vmr_set(vmr_set) == "shared") {
  SHARED_VMR_CANONICAL_POPULATION
} else {
  plot_populations[[1]]
}

base_dir <- here("heritability", "elastic_net_model", cohort,
                 "tissue_comparison", "regulatory_context", "_m")
fig_dir <- regctx_fig_dir(cohort, population, vmr_set)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pal_group <- c(
  "Heritable" = "#3B6EA8",
  "Non-heritable" = "#C65146",
  "Low prediction" = "#7A7A7A"
)

pal_layer_line <- c(
  "Expression · nearest gene" = "#3B6EA8",
  "PSI" = "#5A8F61"
)

tissue_levels <- c("dlpfc", "hippocampus", "caudate")
tissue_labels <- c(
  dlpfc = "DLPFC",
  hippocampus = "Hippocampus",
  caudate = "Caudate"
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
      panel.spacing = unit(0.65, "lines"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

read_many <- function(pattern, populations = plot_populations) {
  files <- list.files(base_dir, pattern = pattern, recursive = TRUE,
                      full.names = TRUE)
  if (length(files) == 0) return(tibble())
  info <- file.info(files)
  files <- files[!is.na(info$size) & info$size > 0]
  if (length(files) == 0) return(tibble())
  tabs <- lapply(files, function(fn) {
    tryCatch(fread(fn), error = function(e) NULL)
  })
  tabs <- Filter(Negate(is.null), tabs)
  if (length(tabs) == 0) return(tibble())
  out <- rbindlist(tabs, fill = TRUE) |> as_tibble()
  if (!"population" %in% names(out)) return(out)
  out |> filter(population %in% populations)
}

## Shared matched VMRs: category-stratified summaries are population-invariant.
dedupe_shared_population_invariant <- function(df, value_cols) {
  if (!nrow(df) || length(plot_populations) <= 1L) return(df)
  key_cols <- setdiff(
    names(df),
    c(value_cols, "population", "cohort", "vmr_set")
  )
  df |>
    group_by(across(all_of(key_cols))) |>
    filter(population == canonical_pop) |>
    ungroup()
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
      modality == "expression_nearest_gene" ~ "Expression · nearest gene",
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

factor_tissue <- function(x) {
  ux <- unique(stats::na.omit(x))
  ord <- c(intersect(tissue_levels, ux), setdiff(ux, tissue_levels))
  factor(x, levels = ord, labels = tissue_labels[ord])
}

main_layers <- function(df) {
  if (!nrow(df) || !"regctx_layer" %in% names(df)) return(df)
  df |> filter(regctx_layer != "Expression · ABC")
}

abc_layers <- function(df) {
  if (!nrow(df) || !"regctx_layer" %in% names(df)) return(tibble())
  df |> filter(regctx_layer == "Expression · ABC")
}

bounded_log2_or <- function(x, cap = 4) {
  y <- log2(x)
  y[is.infinite(y) & y > 0] <- cap
  y[is.infinite(y) & y < 0] <- -cap
  pmax(pmin(y, cap), -cap)
}

env_term_label <- function(env, term) {
  env_map <- c(
    smoking = "Smoking",
    codeine = "Codeine",
    morphine = "Morphine",
    cocaine = "Cocaine",
    ethanol = "Ethanol",
    nicotine = "Nicotine",
    amphetamines = "Amphetamines",
    any_trauma_hx = "Any trauma history",
    education = "Education",
    marital_status = "Marital status"
  )
  base <- ifelse(env %in% names(env_map), env_map[env], gsub("_", " ", env))
  term <- sub("^env_value", "", as.character(term))
  term <- gsub("_", " ", term)
  term <- dplyr::case_when(
    is.na(term) | term == "" | term == "TRUE" ~ "",
    term == "FALSE" ~ "no",
    term == "less than hs" ~ "< high school",
    term == "more than hs" ~ "> high school",
    TRUE ~ term
  )
  ifelse(term == "", base, paste(base, term, sep = ": "))
}

order_env_labels <- function(df, label_col = "env_label") {
  ord <- df |>
    mutate(
      best_fdr = pmin(vmr_fdr, wilcoxon_fdr, fisher_fdr, na.rm = TRUE),
      best_fdr = ifelse(is.infinite(best_fdr), NA_real_, best_fdr),
      evidence_rank = ifelse(is.na(best_fdr), 1, best_fdr),
      shift_rank = ifelse(is.finite(vmr_shift), abs(vmr_shift), 0)
    ) |>
    group_by(.data[[label_col]]) |>
    summarise(
      best_fdr = min(evidence_rank, na.rm = TRUE),
      max_shift = max(shift_rank, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(best_fdr, desc(max_shift), .data[[label_col]]) |>
    pull(.data[[label_col]])
  rev(ord)
}

save_panel <- function(plot, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = 450)
}

group_df <- read_many("h2_group_association_summary.tsv") |>
  ensure_vmr_set_column() |>
  filter(vmr_set == !!vmr_set) |>
  dedupe_shared_population_invariant(
    value_cols = c(
      "n_vmrs", "n_any_fdr05", "n_any_fdr10",
      "prop_any_fdr05", "prop_any_fdr10",
      "median_max_abs_beta", "median_pairs_tested"
    )
  ) |>
  ensure_run_tag() |>
  regctx_plot_layer() |>
  mutate(regctx_layer = factor_layer(regctx_layer))

bin_df <- read_many("h2_bin_association_summary.tsv") |>
  ensure_vmr_set_column() |>
  filter(vmr_set == !!vmr_set, bin_type == "h2_quintile") |>
  ensure_run_tag() |>
  regctx_plot_layer() |>
  mutate(
    regctx_layer = factor_layer(regctx_layer),
    population = factor(population, levels = plot_populations)
  )

env_df <- read_many("environment_convergence_enrichment.tsv") |>
  ensure_vmr_set_column() |>
  filter(vmr_set == !!vmr_set) |>
  dedupe_shared_population_invariant(
    value_cols = c("odds_ratio", "p_value", "fdr", "n_vmrs", "n_sig")
  ) |>
  ensure_run_tag() |>
  regctx_plot_layer() |>
  mutate(regctx_layer = factor_layer(regctx_layer))

env_wilcoxon <- read_many("environment_convergence_wilcoxon.tsv") |>
  ensure_vmr_set_column() |>
  filter(vmr_set == !!vmr_set) |>
  dedupe_shared_population_invariant(
    value_cols = c("location_shift", "p_value", "fdr")
  ) |>
  ensure_run_tag() |>
  regctx_plot_layer() |>
  mutate(regctx_layer = factor_layer(regctx_layer))

env_vmr <- read_many("environment_convergence_vmr_level.tsv") |>
  ensure_vmr_set_column() |>
  filter(vmr_set == !!vmr_set) |>
  dedupe_shared_population_invariant(
    value_cols = c("location_shift", "p_value", "fdr")
  ) |>
  ensure_run_tag() |>
  regctx_plot_layer() |>
  mutate(regctx_layer = factor_layer(regctx_layer))

prox_group <- read_many("proximity_h2_group_summary.tsv") |>
  ensure_vmr_set_column() |>
  filter(vmr_set == !!vmr_set) |>
  dedupe_shared_population_invariant(
    value_cols = c(
      "n_vmrs", "median_distance", "q25_distance", "q75_distance",
      "prop_within_10kb", "prop_within_100kb"
    )
  )
prox_bin <- read_many("proximity_h2_bin_summary.tsv") |>
  ensure_vmr_set_column() |>
  filter(vmr_set == !!vmr_set, bin_type == "h2_quintile") |>
  mutate(population = factor(population, levels = plot_populations))

plots <- list()

group_main <- main_layers(group_df)
bin_main <- main_layers(bin_df)
env_main <- main_layers(env_df)
env_wilcoxon_main <- main_layers(env_wilcoxon)
env_vmr_main <- main_layers(env_vmr)

if (nrow(group_main) > 0) {
  p_group <- group_main |>
    mutate(
      h2_category = factor(h2_category,
                           levels = c("Heritable", "Non-heritable", "Low prediction")),
      tissue = factor_tissue(tissue)
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

if (nrow(bin_main) > 0) {
  p_bin <- bin_main |>
    mutate(tissue = factor_tissue(tissue)) |>
    ggplot(aes(h2_bin, prop_any_fdr10,
               group = interaction(regctx_layer, population),
               color = regctx_layer,
               linetype = population)) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.8) +
    facet_wrap(~ tissue, nrow = 1) +
    scale_color_manual(values = pal_layer_line, drop = FALSE) +
    scale_linetype_manual(values = c(AA = "solid", EA = "22"), drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = "h2 quintile", y = "VMRs with linked feature",
         color = NULL, linetype = "Population") +
    theme_regctx() +
    theme(legend.position = "top")
  plots$bin <- p_bin
  save_panel(p_bin, "panel_h2_quintile_burden", 6.8, 2.8)
}

env_keys <- c(
  "cohort", "tissue", "population", "modality", "run_tag",
  "regctx_layer", "env", "term"
)

env_evidence <- tibble()
if (nrow(env_main) > 0 || nrow(env_wilcoxon_main) > 0 || nrow(env_vmr_main) > 0) {
  fisher_part <- env_main |>
    select(any_of(env_keys), fisher_or = odds_ratio, fisher_fdr = fdr,
           fisher_p = p_value)
  wilcoxon_part <- env_wilcoxon_main |>
    select(any_of(env_keys), wilcoxon_shift = location_shift,
           wilcoxon_fdr = fdr, wilcoxon_p = p_value)
  vmr_part <- env_vmr_main |>
    select(any_of(env_keys), vmr_shift = location_shift,
           vmr_fdr = fdr, vmr_p = p_value)

  env_evidence <- full_join(vmr_part, wilcoxon_part, by = env_keys) |>
    full_join(fisher_part, by = env_keys) |>
    mutate(
      regctx_layer = factor_layer(regctx_layer),
      tissue = factor_tissue(tissue),
      tissue_num = as.integer(tissue),
      env_label = env_term_label(env, term),
      fisher_sig = !is.na(fisher_fdr) & fisher_fdr < 0.10,
      wilcoxon_sig = !is.na(wilcoxon_fdr) & wilcoxon_fdr < 0.10,
      vmr_sig = !is.na(vmr_fdr) & vmr_fdr < 0.10
    )
}

if (nrow(env_evidence) > 0) {
  env_levels <- order_env_labels(env_evidence)
  shift_lim <- max(abs(env_evidence$vmr_shift), na.rm = TRUE)
  if (!is.finite(shift_lim) || shift_lim == 0) shift_lim <- 0.1
  shift_lim <- ceiling(shift_lim * 20) / 20

  evidence_marks <- bind_rows(
    env_evidence |>
      transmute(tissue_num = tissue_num - 0.22, env_label, regctx_layer,
                test = "VMR-level", sig = vmr_sig),
    env_evidence |>
      transmute(tissue_num = tissue_num, env_label, regctx_layer,
                test = "Feature-rank", sig = wilcoxon_sig),
    env_evidence |>
      transmute(tissue_num = tissue_num + 0.22, env_label, regctx_layer,
                test = "Thresholded", sig = fisher_sig)
  ) |>
    filter(sig) |>
    mutate(
      env_label = factor(env_label, levels = env_levels),
      test = factor(test, levels = c("VMR-level", "Feature-rank", "Thresholded"))
    )

  p_env <- env_evidence |>
    mutate(env_label = factor(env_label, levels = env_levels)) |>
    ggplot(aes(tissue_num, env_label)) +
    geom_tile(aes(fill = vmr_shift), color = "white",
              linewidth = 0.3, width = 0.92, height = 0.9) +
    geom_point(data = evidence_marks, aes(tissue_num, env_label, shape = test),
               color = "black", fill = "black", size = 1.65,
               stroke = 0.35, inherit.aes = FALSE) +
    facet_wrap(~ regctx_layer, nrow = 1) +
    scale_x_continuous(breaks = seq_along(tissue_levels),
                       labels = tissue_labels[tissue_levels],
                       expand = expansion(add = 0.05)) +
    scale_fill_gradient2(
      low = "#3B6EA8", mid = "white", high = "#C65146",
      midpoint = 0, limits = c(-shift_lim, shift_lim),
      na.value = "grey92", name = "VMR-level shift"
    ) +
    scale_shape_manual(
      values = c("VMR-level" = 21, "Feature-rank" = 24, "Thresholded" = 22),
      name = "FDR < 0.10"
    ) +
    labs(x = NULL, y = NULL) +
    theme_regctx() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "right",
      panel.grid = element_blank()
    )
  plots$env <- p_env
  save_panel(p_env, "panel_environment_convergence", 6.8, 4.2)

  p_env_fisher <- env_main |>
    mutate(
      env_label = env_term_label(env, term),
      env_label = factor(env_label, levels = env_levels),
      tissue = factor_tissue(tissue),
      log_or = bounded_log2_or(odds_ratio),
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
  save_panel(p_env_fisher, "supplement_environment_thresholded_overlap", 6.8, 3.8)
}

if (nrow(prox_group) > 0) {
  p_prox_group <- prox_group |>
    mutate(
      h2_category = factor(h2_category,
                           levels = c("Heritable", "Non-heritable", "Low prediction")),
      feature_type = factor(feature_type, levels = c("gene", "psi")),
      tissue = factor_tissue(tissue)
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
  save_panel(p_prox_group, "panel_nearest_feature_h2_group", 4.8, 3.4)
}

if (nrow(prox_bin) > 0) {
  p_prox_bin <- prox_bin |>
    mutate(
      feature_type = factor(feature_type, levels = c("gene", "psi")),
      tissue = factor_tissue(tissue)
    ) |>
    ggplot(aes(h2_bin, median_distance + 1,
               group = interaction(feature_type, population),
               color = feature_type,
               linetype = population)) +
    geom_line(linewidth = 0.55) +
    geom_point(size = 1.7) +
    facet_wrap(~ tissue, nrow = 1) +
    scale_color_manual(values = c(gene = "#3B6EA8", psi = "#5A8F61")) +
    scale_linetype_manual(values = c(AA = "solid", EA = "22"), drop = FALSE) +
    scale_y_log10(labels = label_number(scale_cut = cut_short_scale())) +
    labs(x = "h2 quintile", y = "Median nearest distance, bp",
         linetype = "Population") +
    theme_regctx() +
    theme(legend.position = "top")
  plots$prox_bin <- p_prox_bin
  save_panel(p_prox_bin, "panel_nearest_feature_h2_quintile", 6.8, 2.4)
}

abc_group <- abc_layers(group_df)
abc_bin <- abc_layers(bin_df)
abc_env <- abc_layers(env_df)

if (nrow(abc_group) > 0) {
  p_abc_group <- abc_group |>
    mutate(
      h2_category = factor(h2_category,
                           levels = c("Heritable", "Non-heritable", "Low prediction")),
      tissue = factor_tissue(tissue)
    ) |>
    ggplot(aes(h2_category, prop_any_fdr10, fill = h2_category)) +
    geom_col(width = 0.68, color = "white", linewidth = 0.2) +
    facet_wrap(~ tissue, nrow = 1) +
    scale_fill_manual(values = pal_group, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "VMRs with linked feature") +
    theme_regctx() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "none")
  save_panel(p_abc_group, "supplement_abc_h2_group_burden", 6.8, 2.2)
}

if (nrow(abc_bin) > 0) {
  if (length(plot_populations) > 1L) {
    p_abc_bin <- abc_bin |>
      mutate(tissue = factor_tissue(tissue)) |>
      ggplot(aes(h2_bin, prop_any_fdr10,
                 group = interaction(tissue, population),
                 color = tissue,
                 linetype = population)) +
      geom_line(linewidth = 0.55) +
      geom_point(size = 1.8) +
      scale_color_manual(values = c(
        DLPFC = "#3B6EA8",
        Hippocampus = "#5A8F61",
        Caudate = "#C65146"
      ), drop = FALSE) +
      scale_linetype_manual(values = c(AA = "solid", EA = "22"), drop = FALSE) +
      scale_y_continuous(labels = percent_format(accuracy = 1),
                         expand = expansion(mult = c(0.02, 0.08))) +
      labs(x = "h2 quintile", y = "VMRs with ABC-linked gene",
           linetype = "Population") +
      theme_regctx() +
      theme(legend.position = "top")
  } else {
    p_abc_bin <- abc_bin |>
      mutate(tissue = factor_tissue(tissue)) |>
      ggplot(aes(h2_bin, prop_any_fdr10, group = tissue, color = tissue)) +
      geom_line(linewidth = 0.55) +
      geom_point(size = 1.8) +
      scale_color_manual(values = c(
        DLPFC = "#3B6EA8",
        Hippocampus = "#5A8F61",
        Caudate = "#C65146"
      ), drop = FALSE) +
      scale_y_continuous(labels = percent_format(accuracy = 1),
                         expand = expansion(mult = c(0.02, 0.08))) +
      labs(x = "h2 quintile", y = "VMRs with ABC-linked gene") +
      theme_regctx() +
      theme(legend.position = "top")
  }
  save_panel(p_abc_bin, "supplement_abc_h2_quintile_burden", 4.8, 2.6)
}

if (nrow(abc_env) > 0) {
  p_abc_env <- abc_env |>
    mutate(
      env_label = env_term_label(env, term),
      tissue = factor_tissue(tissue),
      log_or = bounded_log2_or(odds_ratio),
      sig = fdr < 0.10
    ) |>
    ggplot(aes(env_label, tissue, fill = log_or, alpha = sig)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#C65146",
                         midpoint = 0, limits = c(-4, 4), na.value = "grey90",
                         name = "log2 OR") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.45), guide = "none") +
    labs(x = NULL, y = NULL) +
    theme_regctx() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "right")
  save_panel(p_abc_env, "supplement_abc_environment_thresholded_overlap", 4.8, 2.2)
}

if (length(plots) >= 2) {
  ordered <- plots[c("group", "bin", "prox_group", "prox_bin", "env")]
  ordered <- ordered[!vapply(ordered, is.null, logical(1))]
  panel_heights <- c(
    group = 3.6,
    bin = 2.4,
    prox_group = 3.0,
    prox_bin = 2.2,
    env = 4.1
  )
  main <- wrap_plots(
    ordered, ncol = 1, guides = "collect",
    heights = unname(panel_heights[names(ordered)])
  ) &
    theme(legend.position = "bottom")
  save_panel(main, "main_regulatory_context", 7.2,
             max(4, sum(panel_heights[names(ordered)]) + 0.5))
}

message("Saved regulatory-context figures to: ", fig_dir)

#### Reproducibility ####
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
