#### Repeat Element Enrichment — Manuscript-Quality Figures ####
##
## Reads summary TSVs from 01.overlap_repeat_elements.R and generates:
##
## Main figure
##   A. Forest-plot matrix of logistic regression ORs for repeat class
##      enrichment with continuous h² (intergenic VMRs, adjusted for
##      VMR length and SNP count)
##   B. h2-quintile trajectory of repeat overlap (intergenic VMRs)
##
## Supplemental figures
##   - Family-level logistic regression forest plot (intergenic VMRs)
##   - Fisher's exact forest plot (unadjusted; sensitivity/transparency)

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

POPULATIONS <- c("AA", "EA")

OUT_DIR <- here(
  "heritability", "elastic_net_model", "all_individuals",
  "tissue_comparison", "annotation", "repeat_elements", "_m"
)

## Display settings

TISSUE_ORDER     <- c("Caudate", "DLPFC", "Hippocampus")
COMPARISON_ORDER <- c("All VMRs", "Intergenic VMRs")

CLASS_ORDER <- c("any_repeat", "LINE", "SINE", "LTR", "DNA",
                 "Satellite", "Simple_repeat", "Low_complexity")
CLASS_LABELS <- c(
  any_repeat     = "Any repeat",
  LINE           = "LINE",
  SINE           = "SINE",
  LTR            = "LTR",
  DNA            = "DNA transposon",
  Satellite      = "Satellite",
  Simple_repeat  = "Simple repeat",
  Low_complexity = "Low complexity"
)

FAMILY_ORDER <- c("L1", "Alu", "MIR", "ERV1", "ERVL", "ERVL_MaLR", "ERVK", "SVA")
FAMILY_LABELS <- c(
  L1         = "L1 (LINE)",
  Alu        = "Alu (SINE)",
  MIR        = "MIR (SINE)",
  ERV1       = "ERV1 (LTR)",
  ERVL       = "ERVL (LTR)",
  ERVL_MaLR  = "ERVL-MaLR (LTR)",
  ERVK       = "ERVK (LTR)",
  SVA        = "SVA (Retroposon)"
)

QUINTILE_ANNOTATIONS <- c("in_any_repeat", "in_LINE", "in_SINE",
                          "in_LTR", "in_DNA")
QUINTILE_LABELS <- c(
  in_any_repeat = "Any repeat",
  in_LINE       = "LINE",
  in_SINE       = "SINE",
  in_LTR        = "LTR",
  in_DNA        = "DNA transposon"
)
QUINTILE_COLORS <- c(
  "Any repeat"     = "#3F3F3F",
  "LINE"           = "#B6523A",
  "SINE"           = "#5B8FA8",
  "LTR"            = "#6FA287",
  "DNA transposon" = "#C17C59"
)

# Logistic regression (continuous h²) — used in main figure and family supplemental
LOGISTIC_EFFECT_COLORS <- c(
  "Enriched with h\u00b2" = "#B6523A",
  "Depleted with h\u00b2"  = "#3A6F8F",
  "Not significant"        = "#B6B6B6"
)

# Fisher's exact (heritable vs non-heritable) — used in Fisher's supplemental
FISHER_EFFECT_COLORS <- c(
  "Higher in heritable"     = "#B6523A",
  "Higher in non-heritable" = "#3A6F8F",
  "Not significant"         = "#B6B6B6"
)

BASE_THEME <- theme_classic(base_size = 10) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 10),
    axis.title         = element_text(size = 9),
    axis.text          = element_text(size = 8, colour = "black"),
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 8),
    legend.key.size    = grid::unit(0.42, "cm"),
    panel.border       = element_blank(),
    plot.margin        = margin(5, 8, 5, 8)
  )

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

effect_direction_logistic <- function(log2_or, fdr) {
  case_when(
    fdr < 0.05 & log2_or > 0 ~ "Enriched with h\u00b2",
    fdr < 0.05 & log2_or < 0 ~ "Depleted with h\u00b2",
    TRUE                     ~ "Not significant"
  )
}

effect_direction_fisher <- function(log2_or, fdr) {
  case_when(
    fdr < 0.05 & log2_or > 0 ~ "Higher in heritable",
    fdr < 0.05 & log2_or < 0 ~ "Higher in non-heritable",
    TRUE                     ~ "Not significant"
  )
}

ac_ci <- function(n_overlap, n_vmrs, z = qnorm(0.975)) {
  n_tilde <- n_vmrs + z^2
  p_tilde <- (n_overlap + z^2 / 2) / n_tilde
  margin  <- z * sqrt(p_tilde * (1 - p_tilde) / n_tilde)
  tibble(
    ci_lo = pmax(0, p_tilde - margin),
    ci_hi = pmin(1, p_tilde + margin)
  )
}

for (pop in POPULATIONS) {

  cat(sprintf("Processing population: %s \n", pop))

  # Set population-specific file names

  fisher_fn   <- file.path(OUT_DIR, paste0("fishers_repeat_enrichment_", pop, ".tsv"))
  logistic_fn <- file.path(OUT_DIR, paste0("logistic_repeat_enrichment_", pop, ".tsv"))
  vmr_fn      <- file.path(OUT_DIR, paste0("vmr_repeat_overlap_", pop, ".tsv"))

  ## ============================================================
  ## Figure A — Logistic regression forest plot (class level, intergenic VMRs)
  ##   OR = effect of continuous h² on repeat overlap, adjusted for
  ##   VMR length and SNP count.
  ## ============================================================

  logistic_raw <- fread(logistic_fn)

  logistic_class <- logistic_raw |>
    dplyr::filter(
      comparison == "Intergenic_VMRs",
      annotation %in% CLASS_ORDER
    ) |>
    mutate(
      tissue = factor(tissue, levels = TISSUE_ORDER),
      annotation = factor(annotation, levels = CLASS_ORDER),
      annotation_label = factor(
        recode(as.character(annotation), !!!CLASS_LABELS),
        levels = rev(unname(CLASS_LABELS[CLASS_ORDER]))
      ),
      log2_or    = log2(or),
      ci_lo_plot = log2(ci_lo_or),
      ci_hi_plot = log2(ci_hi_or),
      effect     = effect_direction_logistic(log2_or, fdr)
    ) |>
    dplyr::filter(!is.na(annotation_label))

  logistic_ci_vals <- c(logistic_class$ci_lo_plot, logistic_class$ci_hi_plot)
  logistic_limit   <- max(abs(logistic_ci_vals[is.finite(logistic_ci_vals)]), na.rm = TRUE)
  logistic_limit   <- ceiling(logistic_limit * 2) / 2

  fig_A <- ggplot(
    logistic_class,
    aes(
      x      = log2_or,
      y      = annotation_label,
      xmin   = ci_lo_plot,
      xmax   = ci_hi_plot,
      colour = effect
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      linewidth  = 0.45,
      colour     = "grey55"
    ) +
    geom_errorbar(
      aes(xmin = ci_lo_plot, xmax = ci_hi_plot),
      orientation = "y",
      width       = 0.18,
      linewidth   = 0.55,
      alpha       = 0.95
    ) +
    geom_point(size = 2.5) +
    facet_wrap(~ tissue, nrow = 1) +
    scale_colour_manual(
      values = LOGISTIC_EFFECT_COLORS,
      breaks = c("Enriched with h\u00b2", "Depleted with h\u00b2", "Not significant"),
      name   = NULL
    ) +
    scale_x_continuous(
      limits = c(0, logistic_limit),
      breaks = pretty(c(0, logistic_limit), n = 5),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(
      x = expression(log[2]~OR~"(per unit h"^2*", adjusted)"),
      y = NULL
    ) +
    BASE_THEME +
    theme(
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.35),
      panel.grid.major.y = element_blank(),
      legend.position    = "top",
      legend.direction   = "horizontal",
      axis.text.y        = element_text(size = 8, lineheight = 0.95)
    )

  ## ============================================================
  ## Figure B — h²-quintile trajectory (intergenic VMRs)
  ## ============================================================

  vmr_overlap <- fread(vmr_fn)

  quintile_summary <- vmr_overlap |>
    filter(hg38_genes_intergenic == 1) |>
    mutate(
      tissue = factor(tissue, levels = TISSUE_ORDER)
    ) |>
    group_by(tissue) |>
    mutate(
      h2_quintile = factor(
        paste0("Q", ntile(h2_unscaled, 5)),
        levels = paste0("Q", 1:5)
      )
    ) |>
    ungroup() |>
    select(tissue, h2_quintile, all_of(QUINTILE_ANNOTATIONS)) |>
    pivot_longer(
      cols      = all_of(QUINTILE_ANNOTATIONS),
      names_to  = "annotation",
      values_to = "overlap"
    ) |>
    group_by(tissue, h2_quintile, annotation) |>
    summarise(
      n_vmrs    = n(),
      n_overlap = sum(overlap, na.rm = TRUE),
      frac      = n_overlap / n_vmrs,
      .groups   = "drop"
    ) |>
    (\(df) bind_cols(df, ac_ci(df$n_overlap, df$n_vmrs)))() |>
    mutate(
      annotation_label = factor(
        recode(annotation, !!!QUINTILE_LABELS),
        levels = names(QUINTILE_COLORS)
      )
    )

  fwrite(
    quintile_summary,
    file.path(OUT_DIR, "repeat_quintile_summary.tsv"),
    sep = "\t"
  )

  fig_B <- ggplot(
    quintile_summary,
    aes(
      x      = h2_quintile,
      y      = frac,
      group  = annotation_label,
      colour = annotation_label,
      fill   = annotation_label
    )
  ) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.10, colour = NA) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 2.0) +
    facet_wrap(~ tissue, nrow = 1) +
    scale_colour_manual(values = QUINTILE_COLORS, name = NULL) +
    scale_fill_manual(values = QUINTILE_COLORS, name = NULL) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.06))
    ) +
    labs(
      x = expression(paste("h"^2, " quintile (intergenic VMRs)")),
      y = "VMRs overlapping repeat element (%)"
    ) +
    BASE_THEME +
    theme(
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      legend.position    = "right",
      legend.direction   = "vertical"
    )

  ## ============================================================
  ## Supplemental A — Logistic regression forest plot (family level)
  ## ============================================================

  logistic_family <- logistic_raw |>
    dplyr::filter(
      comparison == "Intergenic_VMRs",
      annotation %in% FAMILY_ORDER
    ) |>
    mutate(
      tissue = factor(tissue, levels = TISSUE_ORDER),
      annotation = factor(annotation, levels = FAMILY_ORDER),
      annotation_label = factor(
        recode(as.character(annotation), !!!FAMILY_LABELS),
        levels = rev(unname(FAMILY_LABELS[FAMILY_ORDER]))
      ),
      log2_or    = log2(or),
      ci_lo_plot = log2(ci_lo_or),
      ci_hi_plot = log2(ci_hi_or),
      effect     = effect_direction_logistic(log2_or, fdr)
    ) |>
    dplyr::filter(!is.na(annotation_label))

  logistic_family_ci <- c(logistic_family$ci_lo_plot, logistic_family$ci_hi_plot)
  logistic_family_limit <- max(abs(logistic_family_ci[is.finite(logistic_family_ci)]), na.rm = TRUE)
  logistic_family_limit <- ceiling(logistic_family_limit * 2) / 2

  fig_family <- ggplot(
    logistic_family,
    aes(
      x      = log2_or,
      y      = annotation_label,
      xmin   = ci_lo_plot,
      xmax   = ci_hi_plot,
      colour = effect
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      linewidth  = 0.45,
      colour     = "grey55"
    ) +
    geom_errorbar(
      aes(xmin = ci_lo_plot, xmax = ci_hi_plot),
      orientation = "y",
      width       = 0.18,
      linewidth   = 0.55,
      alpha       = 0.95
    ) +
    geom_point(size = 2.5) +
    facet_wrap(~ tissue, nrow = 1) +
    scale_colour_manual(
      values = LOGISTIC_EFFECT_COLORS,
      breaks = c("Enriched with h\u00b2", "Depleted with h\u00b2", "Not significant"),
      name   = NULL
    ) +
    scale_x_continuous(
      limits = c(
        min(0, floor(min(logistic_family_ci[is.finite(logistic_family_ci)]) * 2) / 2),
        logistic_family_limit
      ),
      breaks = pretty(c(0, logistic_family_limit), n = 5),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(
      x = expression(log[2]~OR~"(per unit h"^2*", adjusted)"),
      y = NULL
    ) +
    BASE_THEME +
    theme(
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.35),
      panel.grid.major.y = element_blank(),
      legend.position    = "top",
      legend.direction   = "horizontal",
      axis.text.y        = element_text(size = 8, lineheight = 0.95)
    )

  ## ============================================================
  ## Supplemental B — Fisher's exact forest plot (unadjusted)
  ##   Shown for transparency; the depletion seen here (OR < 1 for
  ##   LINE/LTR in heritable) is a length confound corrected by
  ##   the logistic regression in Figure A.
  ## ============================================================

  fisher_raw <- fread(fisher_fn)

  fisher_class_supp <- fisher_raw |>
    dplyr::filter(annotation %in% CLASS_ORDER) |>
    mutate(
      tissue = factor(tissue, levels = TISSUE_ORDER),
      comparison_label = factor(
        recode(comparison,
              "All_VMRs"        = "All VMRs",
              "Intergenic_VMRs" = "Intergenic VMRs"),
        levels = COMPARISON_ORDER
      ),
      annotation = factor(annotation, levels = CLASS_ORDER),
      annotation_label = factor(
        recode(as.character(annotation), !!!CLASS_LABELS),
        levels = rev(unname(CLASS_LABELS[CLASS_ORDER]))
      ),
      log2_or    = log2(or),
      ci_lo_plot = log2(ci_lo),
      ci_hi_plot = log2(ci_hi),
      effect     = effect_direction_fisher(log2_or, fdr)
    ) |>
    dplyr::filter(!is.na(annotation_label))

  fisher_ci_vals <- c(fisher_class_supp$ci_lo_plot, fisher_class_supp$ci_hi_plot)
  fisher_limit   <- max(abs(fisher_ci_vals[is.finite(fisher_ci_vals)]), na.rm = TRUE)
  fisher_limit   <- ceiling(fisher_limit * 2) / 2

  fig_fisher_supp <- ggplot(
    fisher_class_supp,
    aes(
      x      = log2_or,
      y      = annotation_label,
      xmin   = ci_lo_plot,
      xmax   = ci_hi_plot,
      colour = effect
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype   = "dashed",
      linewidth  = 0.45,
      colour     = "grey55"
    ) +
    geom_errorbar(
      aes(xmin = ci_lo_plot, xmax = ci_hi_plot),
      orientation = "y",
      width       = 0.18,
      linewidth   = 0.55,
      alpha       = 0.95
    ) +
    geom_point(size = 2.5) +
    facet_grid(
      rows   = vars(comparison_label),
      cols   = vars(tissue),
      switch = "y"
    ) +
    scale_colour_manual(
      values = FISHER_EFFECT_COLORS,
      breaks = c("Higher in heritable", "Higher in non-heritable", "Not significant"),
      name   = NULL
    ) +
    scale_x_continuous(
      limits = c(-fisher_limit, fisher_limit),
      breaks = pretty(c(-fisher_limit, fisher_limit), n = 5),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(
      x = expression(log[2]~odds~ratio~"(heritable / non-heritable, unadjusted)"),
      y = NULL
    ) +
    BASE_THEME +
    theme(
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.35),
      panel.grid.major.y = element_blank(),
      legend.position    = "top",
      legend.direction   = "horizontal",
      strip.placement    = "outside",
      strip.text.y.left  = element_text(angle = 0, hjust = 0),
      axis.text.y        = element_text(size = 8, lineheight = 0.95)
    )

  ## ============================================================
  ## Assemble and save
  ## ============================================================

  fig_main <- fig_A / fig_B +
    plot_layout(heights = c(1.0, 1.0)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 11, face = "bold"))

  save_plot(fig_main, paste0("repeat_elements_main_", pop), width = 10.8, height = 6.5)
  save_plot(fig_A, paste0("repeat_elements_logistic_class_", pop), width = 8.6,  height = 4.5)
  save_plot(fig_B, paste0("repeat_elements_quintiles_", pop), width = 7.8,  height = 3.5)
  save_plot(fig_family, paste0("repeat_elements_logistic_family_", pop), width = 8.6, height = 5.0)
  save_plot(fig_fisher_supp, paste0("repeat_elements_fishers_class_supp_", pop), width = 8.6, height = 6.0)

}

cat("Figures written to:", OUT_DIR, "\n")

#### Reproducibility ####
cat("Reproducibility information:\n")
print(Sys.time())
print(proc.time())
options(width = 120)
if (requireNamespace("sessioninfo", quietly = TRUE)) {
  sessioninfo::session_info()
}
