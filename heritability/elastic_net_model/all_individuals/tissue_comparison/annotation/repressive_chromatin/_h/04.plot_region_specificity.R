#### Region-Specificity Figures — Repressive Heritable VMRs ####
##
## Reads output TSVs from 03.region_specificity.R and generates:
##
##   A. Sharing rates by VMR subgroup across tissues (dot plot)
##   B. Repressive enrichment in tissue-specific vs shared heritable VMRs
##   C. h² concordance by repressive status (paired dot plot)
##
## Style follows theme_overlap() from 06.plot_manuscript_overlap.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

## Paths

OUT_DIR <- here(
  "heritability", "elastic_net_model", "BA_only",
  "tissue_comparison", "annotation", "repressive_chromatin", "_m"
)

SHARING_F     <- file.path(OUT_DIR, "region_specificity_sharing_rates.tsv")
FISHER_F      <- file.path(OUT_DIR, "region_specificity_fishers.tsv")
CONCORDANCE_F <- file.path(OUT_DIR, "region_specificity_h2_concordance.tsv")
CAUDATE_F     <- file.path(OUT_DIR, "region_specificity_caudate_enrichment.tsv")

## Display settings

TISSUE_ORDER <- c("Caudate", "DLPFC", "Hippocampus")

SUBGROUP_ORDER <- c(
  "repressive_heritable",
  "nonrepressive_heritable",
  "nonheritable"
)
SUBGROUP_LABELS <- c(
  "repressive_heritable"    = "Repressive heritable",
  "nonrepressive_heritable" = "Non-repressive heritable",
  "nonheritable"            = "Non-heritable"
)
SUBGROUP_COLORS <- c(
  "Repressive heritable"     = "#B6523A",
  "Non-repressive heritable" = "#3A6F8F",
  "Non-heritable"            = "#B6B6B6"
)

SUBGROUP_IG_ORDER <- c(
  "repressive_heritable_intergenic",
  "nonrepressive_heritable_intergenic",
  "nonheritable_intergenic"
)
SUBGROUP_IG_LABELS <- c(
  "repressive_heritable_intergenic"    = "Repressive heritable",
  "nonrepressive_heritable_intergenic" = "Non-repressive heritable",
  "nonheritable_intergenic"            = "Non-heritable"
)

SHARING_LABELS <- c(
  "tissue_specific" = "Tissue-specific",
  "shared_any"      = "Shared (2+)",
  "shared_all3"     = "All 3 tissues"
)
SHARING_COLORS <- c(
  "Tissue-specific" = "#D49A72",
  "Shared (2+)"     = "#2F5964",
  "All 3 tissues"   = "#1A3A42"
)

CONCORDANCE_SUBSETS <- c(
  "all_heritable",
  "repressive_both",
  "nonrepressive_both"
)
CONCORDANCE_LABELS <- c(
  "all_heritable"      = "All heritable",
  "repressive_both"    = "Repressive (both)",
  "nonrepressive_both" = "Non-repressive (both)"
)
CONCORDANCE_COLORS <- c(
  "All heritable"           = "#76888A",
  "Repressive (both)"       = "#B6523A",
  "Non-repressive (both)"   = "#3A6F8F"
)

PAIR_LEVELS <- c("Caudate-DLPFC", "Caudate-Hippocampus", "Hippocampus-DLPFC")

## Theme (matching theme_overlap from manuscript overlap figures)

theme_region <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.title         = element_text(face = "bold"),
      axis.text          = element_text(color = "#1A1A1A"),
      strip.text         = element_text(face = "bold", color = "#1A1A1A"),
      strip.background   = element_blank(),
      plot.title         = element_blank(),
      plot.subtitle      = element_blank(),
      legend.title       = element_text(face = "bold"),
      legend.position    = "right"
    )
}

save_plot <- function(plot_obj, name, width, height) {
  ggsave(
    filename = file.path(OUT_DIR, paste0(name, ".pdf")),
    plot     = plot_obj,
    width    = width,
    height   = height,
    useDingbats = FALSE
  )
  ggsave(
    filename = file.path(OUT_DIR, paste0(name, ".png")),
    plot     = plot_obj,
    width    = width,
    height   = height,
    dpi      = 300
  )
}

## Load data

sharing_df     <- fread(SHARING_F)
fisher_df      <- fread(FISHER_F)
concordance_df <- fread(CONCORDANCE_F)
caudate_df     <- fread(CAUDATE_F)

## ============================================================
## Panel A: Sharing rates by VMR subgroup
## ============================================================

panel_a_data <- sharing_df |>
  filter(subgroup %in% SUBGROUP_ORDER) |>
  mutate(
    tissue = factor(tissue, levels = TISSUE_ORDER),
    subgroup_label = factor(
      recode(subgroup, !!!SUBGROUP_LABELS),
      levels = rev(unname(SUBGROUP_LABELS[SUBGROUP_ORDER]))
    )
  )

panel_a <- ggplot(
  panel_a_data,
  aes(
    x     = sharing_rate,
    y     = subgroup_label,
    xmin  = ci_lo,
    xmax  = ci_hi,
    colour = subgroup_label
  )
) +
  geom_vline(
    xintercept = 0.5,
    linetype   = "dashed",
    linewidth  = 0.4,
    colour     = "grey60"
  ) +
  geom_errorbarh(
    height   = 0.2,
    linewidth = 0.55
  ) +
  geom_point(size = 2.8) +
  geom_text(
    aes(label = sprintf("%.0f%%", 100 * sharing_rate)),
    hjust  = -0.3,
    size   = 2.8,
    show.legend = FALSE
  ) +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = SUBGROUP_COLORS, guide = "none") +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = "VMRs shared with \u22651 other tissue (%)",
    y = NULL
  ) +
  theme_region() +
  theme(
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    legend.position    = "none"
  )

## Intergenic supplemental
panel_a_ig_data <- sharing_df |>
  filter(subgroup %in% SUBGROUP_IG_ORDER) |>
  mutate(
    tissue = factor(tissue, levels = TISSUE_ORDER),
    subgroup_label = factor(
      recode(subgroup, !!!SUBGROUP_IG_LABELS),
      levels = rev(unname(SUBGROUP_IG_LABELS[SUBGROUP_IG_ORDER]))
    )
  )

panel_a_ig <- ggplot(
  panel_a_ig_data,
  aes(
    x     = sharing_rate,
    y     = subgroup_label,
    xmin  = ci_lo,
    xmax  = ci_hi,
    colour = subgroup_label
  )
) +
  geom_vline(
    xintercept = 0.5,
    linetype   = "dashed",
    linewidth  = 0.4,
    colour     = "grey60"
  ) +
  geom_errorbarh(
    height   = 0.2,
    linewidth = 0.55
  ) +
  geom_point(size = 2.8) +
  geom_text(
    aes(label = sprintf("%.0f%%", 100 * sharing_rate)),
    hjust  = -0.3,
    size   = 2.8,
    show.legend = FALSE
  ) +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_manual(values = SUBGROUP_COLORS, guide = "none") +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = "Intergenic VMRs shared with \u22651 other tissue (%)",
    y = NULL
  ) +
  theme_region() +
  theme(
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    legend.position    = "none"
  )

## ============================================================
## Panel B: Repressive enrichment by sharing status
## ============================================================

panel_b_data <- caudate_df |>
  mutate(
    tissue  = factor(tissue, levels = TISSUE_ORDER),
    sharing_label = factor(
      recode(sharing, !!!SHARING_LABELS),
      levels = unname(SHARING_LABELS)
    )
  )

panel_b <- ggplot(
  panel_b_data,
  aes(
    x    = sharing_label,
    y    = frac,
    ymin = ci_lo,
    ymax = ci_hi,
    fill = sharing_label
  )
) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_errorbar(
    width     = 0.2,
    linewidth = 0.5,
    colour    = "#1A1A1A"
  ) +
  geom_text(
    aes(label = sprintf("%.0f%%\nn=%s", 100 * frac, comma(n_total))),
    vjust     = -0.5,
    size      = 2.6,
    lineheight = 0.9
  ) +
  facet_wrap(~ tissue, nrow = 1) +
  scale_fill_manual(values = SHARING_COLORS, guide = "none") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.25))
  ) +
  labs(
    x = "Heritable VMR sharing status",
    y = "Fraction overlapping\nbroad repressive chromatin"
  ) +
  theme_region() +
  theme(
    axis.text.x     = element_text(size = 8, angle = 25, hjust = 1),
    legend.position = "none"
  )

## ============================================================
## Panel C: h² concordance by repressive status
## ============================================================

panel_c_data <- concordance_df |>
  filter(subset %in% CONCORDANCE_SUBSETS) |>
  mutate(
    pair_label = factor(pair_label, levels = PAIR_LEVELS),
    subset_label = factor(
      recode(subset, !!!CONCORDANCE_LABELS),
      levels = unname(CONCORDANCE_LABELS[CONCORDANCE_SUBSETS])
    )
  )

panel_c <- ggplot(
  panel_c_data,
  aes(
    x      = pair_label,
    y      = spearman_rho,
    colour = subset_label,
    group  = subset_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype   = "dashed",
    linewidth  = 0.4,
    colour     = "grey60"
  ) +
  geom_line(
    linewidth  = 0.6,
    alpha      = 0.5
  ) +
  geom_point(
    aes(size = n_pairs),
    alpha = 0.85
  ) +
  geom_text(
    aes(label = sprintf("\u03C1=%.2f\nn=%s", spearman_rho, comma(n_pairs))),
    vjust       = -1.2,
    size        = 2.5,
    show.legend = FALSE,
    lineheight  = 0.9
  ) +
  scale_colour_manual(
    values = CONCORDANCE_COLORS,
    name   = "VMR subset"
  ) +
  scale_size_continuous(
    range = c(2, 5),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(-0.1, 1.0),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(
    x = "Tissue pair",
    y = expression("Spearman " * rho * " (h"^2 * " concordance)")
  ) +
  theme_region() +
  theme(
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    legend.position    = "top",
    legend.direction   = "horizontal"
  )

## ============================================================
## Assemble and save
## ============================================================

fig_main <- (panel_a / panel_b / panel_c) +
  plot_layout(heights = c(0.8, 1.0, 1.0)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold"))

save_plot(fig_main, "region_specificity_main", width = 10.5, height = 10.5)
save_plot(panel_a, "region_specificity_sharing_rates", width = 9.5, height = 3.0)
save_plot(panel_a_ig, "region_specificity_sharing_rates_intergenic", width = 9.5, height = 3.0)
save_plot(panel_b, "region_specificity_repressive_by_sharing", width = 9.5, height = 3.8)
save_plot(panel_c, "region_specificity_h2_concordance", width = 7.5, height = 4.5)

cat("Figures written to:", OUT_DIR, "\n")

#### Reproducibility ####
cat("\nReproducibility information:\n")
print(Sys.time())
print(proc.time())
options(width = 120)
if (requireNamespace("sessioninfo", quietly = TRUE)) {
  sessioninfo::session_info()
}
