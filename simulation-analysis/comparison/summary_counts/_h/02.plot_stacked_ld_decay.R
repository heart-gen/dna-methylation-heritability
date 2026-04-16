# Manuscript-ready phenotype counts across LD decay for two elastic-net methods
# and one target reference, comparing boosting_hybrid vs. joint_ridge portability

suppressPackageStartupMessages({
  library(ggplot2)
  library(here)
  library(dplyr)
  library(tidyr)
  library(readr)
})

## Configuration ##
ld_levels        <- c("0.8", "0.7", "0.6", "0.5")
enet_methods     <- c("boosting_hybrid", "joint_ridge")
method_labels    <- c(boosting_hybrid = "Boosting-hybrid", joint_ridge = "Joint-ridge")
category_levels  <- c("Heritable", "Non-heritable", "Low prediction")
facet_levels     <- c("Target", "Boosting-hybrid", "Joint-ridge")
category_colors  <- c(
  "Heritable"       = "#2B6F8A",
  "Non-heritable"   = "#7E9B63",
  "Low prediction"  = "#D38A5C"
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
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      plot.margin = margin(8, 10, 8, 10)
    )
}

save_plot <- function(p, fn, w, h, dpi = 600) {
  for (ext in c(".png", ".pdf")) {
    ggsave(
      file = paste0(fn, ext),
      plot = p,
      width = w,
      height = h,
      dpi = dpi,
      bg = "white"
    )
  }
}

## Data helpers ##
get_enet_path <- function(ld_decay, enet_method) {
  if (ld_decay == "0.8") {
    ## Baseline uses sim_200_indiv (original pipeline, no METHOD label)
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

read_enet_summary <- function(ld_decay, enet_method) {
  enet_path <- get_enet_path(ld_decay, enet_method)
  assert_files_exist(enet_path)

  enet <- read_tsv(enet_path, show_col_types = FALSE)

  if (!"LD_Decay" %in% names(enet)) {
    enet <- enet |> mutate(LD_Decay = ld_decay)
  }

  enet |>
    mutate(
      LD_Decay     = factor(as.character(LD_Decay), levels = ld_levels, ordered = TRUE),
      enet_method  = enet_method,
      facet_method = method_labels[enet_method]
    )
}

read_target_summary <- function() {
  target_path <- get_target_path("0.8")
  assert_files_exist(target_path)
  read_tsv(target_path, show_col_types = FALSE) |>
    mutate(x_group = "Target", facet_method = "Target")
}

## Main ##
out_path <- here("simulation-analysis/comparison/summary_counts/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

all_enet_paths <- unlist(lapply(enet_methods, function(m) {
  vapply(ld_levels, function(ld) get_enet_path(ld, m), FUN.VALUE = character(1))
}))
assert_files_exist(c(all_enet_paths, get_target_path("0.8")))

## Elastic-net summaries for both methods
enet_summary <- bind_rows(lapply(enet_methods, function(m) {
  bind_rows(lapply(ld_levels, function(ld) read_enet_summary(ld, m)))
})) |>
  mutate(
    x_group = as.character(LD_Decay),
    h2_category = case_when(
      !is.finite(r_squared_cv) | !is.finite(h2_unscaled) ~ "Low prediction",
      r_squared_cv <= 0.75                                ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.75            ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.75           ~ "Heritable"
    ),
    h2_category  = factor(h2_category, levels = category_levels),
    facet_method = factor(facet_method, levels = facet_levels)
  ) |>
  count(x_group, facet_method, h2_category, name = "count")

target_summary <- read_target_summary() |>
  mutate(
    h2_category = case_when(
      target_heritability < 0.1  ~ "Non-heritable",
      target_heritability >= 0.1 ~ "Heritable"
    ),
    h2_category  = factor(h2_category, levels = category_levels),
    facet_method = factor(facet_method, levels = facet_levels)
  ) |>
  count(x_group, facet_method, h2_category, name = "count")

x_levels <- c("Target", ld_levels)

all_summary <- bind_rows(enet_summary, target_summary) |>
  mutate(
    facet_method = factor(facet_method, levels = facet_levels),
    x_group      = factor(x_group, levels = x_levels),
    h2_category  = factor(h2_category, levels = category_levels)
  ) |>
  group_by(facet_method, x_group) |>
  complete(h2_category, fill = list(count = 0)) |>
  ungroup() |>
  arrange(facet_method, x_group, h2_category)

count_check <- all_summary |>
  group_by(facet_method, x_group) |>
  summarise(total = sum(count), .groups = "drop")

if (any(count_check$total[!is.na(count_check$total)] != 1000, na.rm = TRUE)) {
  stop("Category counts did not sum to 1000 for all plotted method-LD groups.")
}

print(all_summary)

p <- ggplot(all_summary, aes(x = x_group, y = count, fill = h2_category)) +
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.35
  ) +
  facet_grid(. ~ facet_method, switch = "x", scales = "free_x", space = "free_x") +
  scale_fill_manual(values = category_colors, drop = FALSE) +
  scale_y_continuous(
    breaks = seq(0, 1000, by = 250),
    expand = expansion(mult = c(0, 0.02))
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(
    x = NULL,
    y = "Phenotypes (n)"
  ) +
  publication_theme()

plot_file <- file.path(out_path, "simulated_data_stacked_ld_decay_200")
save_plot(p, plot_file, w = 8.5, h = 4.8)

## Reproducibility ##
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
