#### Publication figures for environmental enrichment discovery and replication ####

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(stringr)
})

find_repo_root <- function(start_dir = getwd()) {
  cur <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(cur, ".git")) ||
        dir.exists(file.path(cur, "environmental-analysis"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      stop("Could not locate repository root from: ", start_dir)
    }
    cur <- parent
  }
}

repo_root <- find_repo_root()
discovery_dir <- file.path(
  repo_root, "environmental-analysis", "BA_only", "tissue_compare",
  "correlation", "enrichment", "_m"
)
replication_dir <- file.path(
  repo_root, "environmental-analysis", "all_individuals", "tissue_compare",
  "correlation", "enrichment", "_m"
)
out_dir <- file.path(discovery_dir, "publication_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cohort_meta <- tibble::tribble(
  ~cohort_key, ~Cohort,             ~cohort_short, ~file,
  "discovery", "Discovery",        "Discovery",   file.path(discovery_dir, "vmr_enrichment_analysis.txt"),
  "all",       "Replication (Shared)", "Shared",   file.path(replication_dir, "vmr_enrichment_analysis-all.txt"),
  "BA",        "Replication (BA)",  "BA",          file.path(replication_dir, "vmr_enrichment_analysis-BA.txt"),
  "WA",        "Replication (WA)",  "WA",          file.path(replication_dir, "vmr_enrichment_analysis-WA.txt")
)

cohort_levels <- cohort_meta$Cohort
cohort_short_levels <- cohort_meta$cohort_short
tissue_levels <- c("Caudate", "DLPFC", "Hippocampus")
test_levels <- c("Logit", "DMR", "Both")
test_labels <- c(Logit = "VMR", DMR = "DMR", Both = "VMR + DMR")
h2_levels <- c("Heritable", "Non-heritable")
h2_labels <- c(Heritable = "Heritable", "Non-heritable" = "Non-heritable")
h2_short_labels <- c(Heritable = "H", "Non-heritable" = "N")
h2_palette <- c("Heritable" = "#497C8A", "Non-heritable" = "#8CA77B")

env_label_map <- c(
  less_than_hs = "< high school",
  more_than_hs = "> high school",
  previously_married = "Previously married"
)

env_short_label_map <- c(
  amphetamines = "Amphet.",
  less_than_hs = "<HS",
  more_than_hs = ">HS",
  previously_married = "Prev. married"
)

label_environment <- function(x) {
  out <- ifelse(
    x %in% names(env_label_map),
    unname(env_label_map[x]),
    stringr::str_to_title(gsub("_", " ", x))
  )
  out
}

label_environment_short <- function(x) {
  out <- ifelse(
    x %in% names(env_short_label_map),
    unname(env_short_label_map[x]),
    stringr::str_to_title(gsub("_", " ", x))
  )
  out
}

label_environment_main_axis <- function(x) {
  recode(
    x,
    "Amphetamines" = "Amphet.",
    "Previously married" = "Prev. married",
    .default = x
  )
}

safe_filename <- function(x) {
  x |>
    tolower() |>
    gsub("[^a-z0-9]+", "_", x = _) |>
    gsub("^_|_$", "", x = _)
}

save_figure <- function(plot, stem, width, height) {
  ggsave(
    file.path(out_dir, paste0(stem, ".pdf")),
    plot = plot, width = width, height = height, units = "in",
    bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0(stem, ".png")),
    plot = plot, width = width, height = height, units = "in",
    dpi = 320, bg = "white"
  )
}

theme_enrichment <- function(base_size = 8) {
  ggpubr::theme_pubr(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = base_size, face = "bold"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.45, "cm"),
      panel.spacing = unit(0.45, "lines"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

read_enrichment <- function(path, cohort_key, cohort_label, cohort_short) {
  if (!file.exists(path)) {
    stop("Missing enrichment input: ", path)
  }
  dt <- data.table::fread(path)
  required <- c("Tissue", "h2_Category", "OR", "PValue", "FDR", "Test", "Env")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  }
  dt |>
    as_tibble() |>
    mutate(
      cohort_key = cohort_key,
      Cohort = cohort_label,
      cohort_short = cohort_short,
      Tissue = factor(Tissue, levels = tissue_levels),
      Test = factor(Test, levels = test_levels),
      test_label = factor(test_labels[as.character(Test)],
                          levels = unname(test_labels[test_levels])),
      Env = as.character(Env),
      env_label = label_environment(Env),
      h2_Category = as.character(h2_Category),
      OR = as.numeric(OR),
      PValue = as.numeric(PValue),
      FDR = as.numeric(FDR),
      log2_or = ifelse(is.finite(OR) & OR > 0, log2(OR), NA_real_),
      h2_key = recode(h2_Category, "Non-heritable" = "Nonheritable",
                      .default = h2_Category)
    )
}

enrichment <- purrr::pmap_dfr(
  cohort_meta,
  \(cohort_key, Cohort, cohort_short, file) {
    read_enrichment(file, cohort_key, Cohort, cohort_short)
  }
)

plot_df <- enrichment |>
  filter(h2_Category %in% h2_levels) |>
  mutate(
    Cohort = factor(Cohort, levels = cohort_levels),
    cohort_short = factor(cohort_short, levels = cohort_short_levels),
    h2_Category = factor(h2_Category, levels = h2_levels),
    h2_label = factor(h2_labels[as.character(h2_Category)],
                      levels = unname(h2_labels[h2_levels])),
    h2_short = factor(h2_short_labels[as.character(h2_Category)],
                      levels = unname(h2_short_labels[h2_levels])),
    sig = !is.na(FDR) & FDR < 0.05
  )

contrast_df <- plot_df |>
  select(
    cohort_key, Cohort, cohort_short, Tissue, Test, test_label, Env, env_label,
    h2_key, OR, FDR, log2_or
  ) |>
  tidyr::pivot_wider(
    names_from = h2_key,
    values_from = c(OR, FDR, log2_or),
    names_glue = "{.value}_{h2_key}"
  ) |>
  mutate(
    contrast_log2or = log2_or_Heritable - log2_or_Nonheritable,
    sig_both = !is.na(FDR_Heritable) & !is.na(FDR_Nonheritable) &
      FDR_Heritable < 0.05 & FDR_Nonheritable < 0.05,
    sig_any = !is.na(FDR_Heritable) & FDR_Heritable < 0.05 |
      !is.na(FDR_Nonheritable) & FDR_Nonheritable < 0.05,
    sig_class = case_when(
      sig_both ~ "Both classes",
      sig_any ~ "One class",
      TRUE ~ "Not significant"
    )
  )

env_ord <- c(
    "Smoking",
    "Nicotine",
    "Cocaine",
    "Morphine",
    "Codeine",
    "Amphetamines",
    "Ethanol",
    "< high school",
    "> high school",
    "Previously married",
    "Single"
  )

if (length(env_order) == 0) {
  env_order <- sort(unique(plot_df$env_label))
}

plot_df <- plot_df |>
  mutate(env_label = factor(env_label, levels = rev(env_order)))

contrast_df <- contrast_df |>
  mutate(env_label = factor(env_label, levels = rev(env_order)))

contrast_lim <- max(abs(contrast_df$contrast_log2or), na.rm = TRUE)
contrast_lim <- ifelse(is.finite(contrast_lim), ceiling(contrast_lim * 2) / 2, 1)
contrast_breaks <- c(-contrast_lim, 0, contrast_lim)
contrast_labels <- c("Non-heritable\nstronger", "No\ndifference", "Heritable\nstronger")

or_lim <- max(abs(plot_df$log2_or), na.rm = TRUE)
or_lim <- ifelse(is.finite(or_lim), ceiling(or_lim * 2) / 2, 1)

cohort_test_levels <- tidyr::expand_grid(
  cohort_short = cohort_short_levels,
  test_label = unname(test_labels[test_levels])
) |>
  mutate(cohort_test = paste(cohort_short, test_label, sep = "\n")) |>
  pull(cohort_test)

heatmap_df <- contrast_df |>
  mutate(
    cohort_test = factor(
      paste(as.character(cohort_short), as.character(test_label), sep = "\n"),
      levels = cohort_test_levels
    ),
    sig_class = factor(sig_class, levels = c("One class", "Both classes"))
  )

p_heatmap_all <- ggplot(heatmap_df, aes(cohort_test, env_label, fill = contrast_log2or)) +
  geom_tile(color = "white", linewidth = 0.25, width = 0.95, height = 0.95) +
  geom_point(
    data = heatmap_df |> filter(sig_any),
    aes(shape = sig_class),
    color = "black", fill = "black", size = 1.35, stroke = 0.35
  ) +
  facet_grid(. ~ Tissue) +
  scale_fill_gradient2(
    low = "#C65146", mid = "white", high = "#3B6EA8",
    midpoint = 0, limits = c(-contrast_lim, contrast_lim),
    breaks = contrast_breaks,
    labels = contrast_labels,
    na.value = "grey90",
    name = "Stronger\nenrichment"
  ) +
  scale_shape_manual(
    values = c("One class" = 1, "Both classes" = 21),
    drop = TRUE,
    name = "FDR < 0.05"
  ) +
  labs(x = NULL, y = NULL) +
  theme_enrichment(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    panel.grid = element_blank()
  )

heatmap_both_df <- contrast_df |>
  filter(as.character(Test) == "Both") |>
  mutate(
    cohort_short = factor(cohort_short, levels = cohort_short_levels),
    sig_class = factor(sig_class, levels = c("One class", "Both classes"))
  )

p_heatmap_both <- ggplot(heatmap_both_df, aes(cohort_short, env_label, fill = contrast_log2or)) +
  geom_tile(color = "white", linewidth = 0.25, width = 0.92, height = 0.95) +
  geom_point(
    data = heatmap_both_df |> filter(sig_any),
    aes(shape = sig_class),
    color = "black", fill = "black", size = 1.45, stroke = 0.35
  ) +
  facet_grid(. ~ Tissue) +
  scale_fill_gradient2(
    low = "#C65146", mid = "white", high = "#3B6EA8",
    midpoint = 0, limits = c(-contrast_lim, contrast_lim),
    breaks = contrast_breaks,
    labels = contrast_labels,
    na.value = "grey90",
    name = "Stronger\nenrichment"
  ) +
  scale_shape_manual(
    values = c("One class" = 1, "Both classes" = 21),
    drop = TRUE,
    name = "FDR < 0.05"
  ) +
  scale_y_discrete(labels = label_environment_main_axis) +
  labs(x = NULL, y = NULL) +
  theme_enrichment(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(4, 4, 4, 14)
  )

make_top_contrasts <- function(test_filter = NULL, include_test = TRUE, n = 12) {
  out <- contrast_df |>
    filter(cohort_key == "discovery", is.finite(contrast_log2or))
  if (!is.null(test_filter)) {
    out <- out |> filter(as.character(Test) %in% test_filter)
  }
  out |>
    arrange(desc(abs(contrast_log2or)), env_label, Tissue, test_label) |>
    slice_head(n = n) |>
    mutate(
      contrast_rank = row_number(),
      tissue_short = recode(
        as.character(Tissue),
        Hippocampus = "Hipp.",
        .default = as.character(Tissue)
      ),
      test_short = recode(
        as.character(test_label),
        "VMR + DMR" = "VMR+DMR",
        .default = as.character(test_label)
      ),
      contrast_label = if (include_test) {
        paste0(label_environment_short(Env), " | ", tissue_short, " | ", test_short)
      } else {
        paste0(label_environment_short(Env), " | ", tissue_short)
      }
    ) |>
    select(Tissue, Test, Env, contrast_rank, contrast_label)
}

make_paired_plot <- function(top_contrasts, base_size = 7) {
  paired_df <- plot_df |>
    semi_join(top_contrasts, by = c("Tissue", "Test", "Env")) |>
    left_join(top_contrasts, by = c("Tissue", "Test", "Env")) |>
    mutate(
      contrast_label = factor(
        contrast_label,
        levels = rev(top_contrasts$contrast_label)
      )
    )

  ggplot(
    paired_df,
    aes(x = log2_or, y = contrast_label)
  ) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_line(
      aes(group = interaction(Cohort, Tissue, Test, Env)),
      color = "grey72", linewidth = 0.35, na.rm = TRUE
    ) +
    geom_point(
      aes(color = h2_Category, shape = sig),
      size = 1.75, stroke = 0.45, na.rm = TRUE
    ) +
    facet_grid(. ~ Cohort) +
    scale_color_manual(values = h2_palette, drop = FALSE, name = NULL) +
    scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 1),
      labels = c("TRUE" = "FDR < 0.05", "FALSE" = "n.s."),
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(-or_lim, or_lim),
      breaks = pretty_breaks(n = 5)
    ) +
    labs(x = "log2(OR)", y = NULL) +
    theme_enrichment(base_size = base_size) +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 4, 4, 14)
    )
}

top_contrasts_all <- make_top_contrasts()
top_contrasts_both <- make_top_contrasts(test_filter = "Both", include_test = FALSE)

p_paired_all <- make_paired_plot(top_contrasts_all, base_size = 7)
p_paired_both <- make_paired_plot(top_contrasts_both, base_size = 8)

main_figure <- ggpubr::ggarrange(
  p_heatmap_both, p_paired_both,
  ncol = 1,
  heights = c(1.2, 1),
  labels = c("A", "B"),
  font.label = list(size = 10, face = "bold")
)

save_figure(
  main_figure,
  "main_environmental_enrichment_discovery_replication",
  width = 12.5,
  height = 10.2
)

supplement_all_tests <- ggpubr::ggarrange(
  p_heatmap_all, p_paired_all,
  ncol = 1,
  heights = c(1.2, 1),
  labels = c("A", "B"),
  font.label = list(size = 10, face = "bold")
)

save_figure(
  supplement_all_tests,
  "supplement_environmental_enrichment_discovery_replication_all_tests",
  width = 12.5,
  height = 10.2
)

atlas_x_levels <- tidyr::expand_grid(
  cohort_short = cohort_short_levels,
  h2_short = unname(h2_short_labels[h2_levels])
) |>
  mutate(cohort_h2 = paste(cohort_short, h2_short, sep = "\n")) |>
  pull(cohort_h2)

atlas_df <- plot_df |>
  mutate(
    cohort_h2 = factor(
      paste(as.character(cohort_short), as.character(h2_short), sep = "\n"),
      levels = atlas_x_levels
    )
  )

p_atlas <- ggplot(atlas_df, aes(cohort_h2, test_label, fill = log2_or)) +
  geom_tile(color = "white", linewidth = 0.2, width = 0.95, height = 0.9) +
  geom_point(
    data = atlas_df |> filter(sig),
    shape = 21, color = "black", fill = "black", size = 0.75, stroke = 0.2
  ) +
  facet_grid(env_label ~ Tissue) +
  scale_fill_gradient2(
    low = "#3B6EA8", mid = "white", high = "#C65146",
    midpoint = 0, limits = c(-or_lim, or_lim), na.value = "grey90",
    name = "log2(OR)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_enrichment(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )

save_figure(
  p_atlas,
  "supplement_environment_enrichment_atlas_discovery_replication",
  width = 13.5,
  height = 16
)

plot_environment_panel <- function(env_id) {
  env_grid <- tidyr::expand_grid(
    Cohort = factor(cohort_levels, levels = cohort_levels),
    Tissue = factor(tissue_levels, levels = tissue_levels),
    Test = factor(test_levels, levels = test_levels),
    h2_Category = factor(h2_levels, levels = h2_levels)
  ) |>
    mutate(
      Env = env_id,
      env_label = label_environment(env_id),
      test_label = factor(test_labels[as.character(Test)],
                          levels = unname(test_labels[test_levels])),
      h2_label = factor(h2_labels[as.character(h2_Category)],
                        levels = unname(h2_labels[h2_levels]))
    )

  env_df <- plot_df |>
    filter(Env == env_id) |>
    select(Cohort, Tissue, Test, h2_Category, log2_or, FDR, sig)

  panel_df <- env_grid |>
    left_join(env_df, by = c("Cohort", "Tissue", "Test", "h2_Category")) |>
    mutate(env_label = factor(env_label, levels = label_environment(env_id)))

  ggplot(panel_df, aes(h2_label, test_label, fill = log2_or)) +
    geom_tile(color = "white", linewidth = 0.25, width = 0.95, height = 0.9) +
    geom_point(
      data = panel_df |> filter(!is.na(sig) & sig),
      shape = 21, color = "black", fill = "black", size = 1, stroke = 0.25
    ) +
    facet_grid(Cohort + env_label ~ Tissue) +
    scale_fill_gradient2(
      low = "#3B6EA8", mid = "white", high = "#C65146",
      midpoint = 0, limits = c(-or_lim, or_lim), na.value = "grey90",
      name = "log2(OR)"
    ) +
    labs(x = NULL, y = NULL) +
    theme_enrichment(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "right",
      panel.grid = element_blank(),
      strip.text.y = element_text(angle = 0, hjust = 0)
    )
}

env_ids <- contrast_df |>
  distinct(Env, env_label) |>
  mutate(env_label_chr = as.character(env_label)) |>
  arrange(factor(env_label_chr, levels = rev(env_order))) |>
  pull(Env)

for (env_id in env_ids) {
  env_plot <- plot_environment_panel(env_id)
  save_figure(
    env_plot,
    paste0("supplement_environment_", safe_filename(env_id),
           "_enrichment_discovery_replication"),
    width = 8.8,
    height = 7.2
  )
}

plot_data_file <- file.path(out_dir, "environmental_enrichment_publication_plot_data.tsv")
data.table::fwrite(
  plot_df |>
    mutate(
      Cohort = as.character(Cohort),
      Tissue = as.character(Tissue),
      Test = as.character(Test),
      h2_Category = as.character(h2_Category),
      env_label = as.character(env_label)
    ) |>
    select(
      Cohort, cohort_key, Tissue, Test, test_label, Env, env_label,
      h2_Category, OR, PValue, FDR, log2_or, sig
    ),
  plot_data_file,
  sep = "\t"
)

message("Saved publication enrichment figures to: ", out_dir)
message("Saved plotting data to: ", plot_data_file)

Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
