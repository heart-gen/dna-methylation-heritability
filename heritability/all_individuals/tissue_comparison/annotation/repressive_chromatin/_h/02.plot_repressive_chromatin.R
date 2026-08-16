#### Repressive Chromatin Enrichment — Manuscript-Quality Figures ####
##
## Reads summary TSVs from 01.overlap_repressive_chromatin.R and generates:
##
## Main figure
##   A. Forest-plot matrix of Fisher enrichment (heritable vs non-heritable)
##   B. h2-quintile trajectory of repressive-chromatin overlap
##
## Supplemental figures
##   - Full enrichment forest plot
##   - Full ChromHMM absolute overlap heatmap + fold-change heatmap
##   - Peak-call sensitivity forest plot for broad vs gapped H3K27me3 / H3K9me3

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
  "tissue_comparison", "annotation", "repressive_chromatin", "_m"
)

## Display settings

TISSUE_ORDER <- c("Caudate", "DLPFC", "Hippocampus")
COMPARISON_ORDER <- c("All VMRs", "Intergenic VMRs")

ANNOTATION_ORDER <- c(
  "Bivalent", "H3K27me3", "Polycomb",
  "BroadRepressive", "H3K9me3", "Het", "Quies"
)
ANNOTATION_LABELS <- c(
  "Bivalent"         = "Bivalent",
  "H3K27me3"         = "H3K27me3",
  "Polycomb"         = "Polycomb\n(ReprPC)",
  "BroadRepressive"  = "Broad repressive",
  "H3K9me3"          = "H3K9me3",
  "Het"              = "Heterochromatin\n(Het)",
  "Quies"            = "Quiescent\n(Quies)"
)

PEAK_TYPE_ORDER <- c("gappedPeak", "broadPeak")
MARK_ORDER <- c("H3K27me3", "H3K9me3")

STATE_ORDER_FULL <- c(
  "1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk",
  "6_EnhG", "7_Enh", "8_ZNF/Rpts",
  "10_TssBiv", "11_BivFlnk", "12_EnhBiv",
  "9_Het", "13_ReprPC", "14_ReprPCWk", "15_Quies"
)
STATE_LABELS <- c(
  "1_TssA"      = "TssA",
  "2_TssAFlnk"  = "TssAFlnk",
  "3_TxFlnk"    = "TxFlnk",
  "4_Tx"        = "Tx",
  "5_TxWk"      = "TxWk",
  "6_EnhG"      = "EnhG",
  "7_Enh"       = "Enh",
  "8_ZNF/Rpts"  = "ZNF/Rpts",
  "10_TssBiv"   = "TssBiv",
  "11_BivFlnk"  = "BivFlnk",
  "12_EnhBiv"   = "EnhBiv",
  "9_Het"       = "Het",
  "13_ReprPC"   = "ReprPC",
  "14_ReprPCWk" = "ReprPCWk",
  "15_Quies"    = "Quies"
)
QUINTILE_ANNOTATIONS <- c(
  "in_BroadRepressive", "in_Quies", "in_H3K9me3",
  "in_H3K27me3", "in_Polycomb"
)
QUINTILE_LABELS <- c(
  "in_BroadRepressive" = "Broad repressive",
  "in_Quies"           = "Quiescent",
  "in_H3K9me3"         = "H3K9me3",
  "in_H3K27me3"        = "H3K27me3",
  "in_Polycomb"        = "Polycomb"
)
QUINTILE_COLORS <- c(
  "Broad repressive" = "#3F3F3F",
  "Quiescent"        = "#C17C59",
  "H3K9me3"          = "#B6523A",
  "H3K27me3"         = "#5B8FA8",
  "Polycomb"         = "#6FA287"
)

EFFECT_COLORS <- c(
  "Higher in heritable"     = "#B6523A",
  "Higher in non-heritable" = "#3A6F8F",
  "Not significant"         = "#B6B6B6"
)
ABS_RATE_COLORS <- c("#F7F3EE", "#D49A72", "#8E4426")

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

effect_direction <- function(log2_or, fdr) {
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

  fisher_fn <- file.path(OUT_DIR, paste0("fishers_repressive_enrichment_", pop, ".tsv"))
  rate_fn   <- file.path(OUT_DIR, paste0("chromhmm_state_overlap_rates_", pop, ".tsv"))
  vmr_fn    <- file.path(OUT_DIR, paste0("vmr_repressive_overlap_", pop, ".tsv"))

  ## Load and prepare Fisher enrichment data

  fisher_raw <- fread(fisher_fn)

  fisher_main <- fisher_raw |>
    filter(!grepl("_broad$", annotation)) |>
    mutate(
      tissue = factor(tissue, levels = TISSUE_ORDER),
      comparison_label = factor(
        recode(comparison,
              "All_VMRs" = "All VMRs",
              "Intergenic_VMRs" = "Intergenic VMRs"),
        levels = COMPARISON_ORDER
      ),
      annotation = factor(annotation, levels = ANNOTATION_ORDER),
      annotation_label = factor(
        recode(as.character(annotation), !!!ANNOTATION_LABELS),
        levels = rev(unname(ANNOTATION_LABELS[ANNOTATION_ORDER]))
      ),
      log2_or = log2(or),
      ci_lo_plot = log2(ci_lo),
      ci_hi_plot = log2(ci_hi),
      effect = effect_direction(log2_or, fdr)
    ) |>
    filter(!is.na(annotation_label))

  fisher_ci_vals <- c(fisher_main$ci_lo_plot, fisher_main$ci_hi_plot)
  fisher_limit <- max(abs(fisher_ci_vals[is.finite(fisher_ci_vals)]), na.rm = TRUE)
  fisher_limit <- ceiling(fisher_limit * 2) / 2

  fig_A <- ggplot(
    fisher_main,
    aes(
      x = log2_or,
      y = annotation_label,
      xmin = ci_lo_plot,
      xmax = ci_hi_plot,
      colour = effect
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.45,
      colour = "grey55"
    ) +
    geom_errorbar(
      aes(xmin = ci_lo_plot, xmax = ci_hi_plot),
      orientation = "y",
      width = 0.18,
      linewidth = 0.55,
      alpha = 0.95
    ) +
    geom_point(size = 2.5) +
    facet_grid(
      rows = vars(comparison_label),
      cols = vars(tissue),
      switch = "y"
    ) +
    scale_colour_manual(
      values = EFFECT_COLORS,
      breaks = c(
        "Higher in heritable",
        "Higher in non-heritable",
        "Not significant"
      ),
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(-fisher_limit, fisher_limit),
      breaks = pretty(c(-fisher_limit, fisher_limit), n = 5),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(
      x = expression(log[2]~odds~ratio~"(heritable / non-heritable)"),
      y = NULL
    ) +
    BASE_THEME +
    theme(
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.35),
      panel.grid.major.y = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 0),
      axis.text.y = element_text(size = 8, lineheight = 0.95)
    )

  ## Sensitivity: peak caller comparison

  fisher_sens <- fisher_raw |>
    filter(
      annotation %in% c("H3K27me3", "H3K9me3", "H3K27me3_broad", "H3K9me3_broad"),
      comparison == "Intergenic_VMRs"
    ) |>
    mutate(
      tissue = factor(tissue, levels = TISSUE_ORDER),
      mark = factor(gsub("_broad$", "", annotation), levels = MARK_ORDER),
      mark_label = factor(mark, levels = rev(MARK_ORDER)),
      peak_type = factor(
        ifelse(grepl("_broad$", annotation), "broadPeak", "gappedPeak"),
        levels = rev(PEAK_TYPE_ORDER)
      ),
      log2_or = log2(or),
      ci_lo_plot = log2(ci_lo),
      ci_hi_plot = log2(ci_hi),
      effect = effect_direction(log2_or, fdr)
    )

  fig_A_sens <- ggplot(
    fisher_sens,
    aes(
      x = log2_or,
      y = peak_type,
      xmin = ci_lo_plot,
      xmax = ci_hi_plot,
      colour = effect
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.45,
      colour = "grey55"
    ) +
    geom_errorbar(
      aes(xmin = ci_lo_plot, xmax = ci_hi_plot),
      orientation = "y",
      width = 0.18,
      linewidth = 0.55
    ) +
    geom_point(size = 2.5) +
    facet_grid(
      rows = vars(mark_label),
      cols = vars(tissue),
      switch = "y"
    ) +
    scale_colour_manual(values = EFFECT_COLORS, guide = "none") +
    scale_x_continuous(
      limits = c(-fisher_limit, fisher_limit),
      breaks = pretty(c(-fisher_limit, fisher_limit), n = 5),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_y_discrete(labels = c("broadPeak" = "broadPeak", "gappedPeak" = "gappedPeak")) +
    labs(
      x = expression(log[2]~odds~ratio~"(heritable / non-heritable)"),
      y = NULL
    ) +
    BASE_THEME +
    theme(
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.35),
      panel.grid.major.y = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 0),
      axis.text.y = element_text(size = 8)
    )

  ## Quintile trajectory summary from per-VMR overlap data

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
      cols = all_of(QUINTILE_ANNOTATIONS),
      names_to = "annotation",
      values_to = "overlap"
    ) |>
    group_by(tissue, h2_quintile, annotation) |>
    summarise(
      n_vmrs = n(),
      n_overlap = sum(overlap, na.rm = TRUE),
      frac = n_overlap / n_vmrs,
      .groups = "drop"
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
    file.path(OUT_DIR, paste0("repressive_quintile_summary_", pop, ".tsv")),
    sep = "\t"
  )

  fig_B_quintile <- ggplot(
    quintile_summary,
    aes(
      x = h2_quintile,
      y = frac,
      group = annotation_label,
      colour = annotation_label,
      fill = annotation_label
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
      y = "VMRs overlapping repressive chromatin (%)"
    ) +
    BASE_THEME +
    theme(
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      legend.position = "right",
      legend.direction = "vertical"
    )

  ## Load and prepare ChromHMM overlap data

  rate_raw <- fread(rate_fn) |>
    mutate(
      tissue = factor(tissue, levels = TISSUE_ORDER),
      h2_category = factor(h2_category, levels = c("Heritable", "Non-heritable")),
      state = factor(state, levels = STATE_ORDER_FULL)
    ) |>
    filter(!is.na(state))

  rate_abs <- rate_raw |>
    mutate(
      state_label = factor(
        recode(as.character(state), !!!STATE_LABELS),
        levels = rev(unname(STATE_LABELS[STATE_ORDER_FULL]))
      )
    )

  rate_diff <- rate_raw |>
    select(tissue, state, h2_category, frac) |>
    pivot_wider(names_from = h2_category, values_from = frac) |>
    mutate(
      log2fc = log2((Heritable + 1e-4) / (`Non-heritable` + 1e-4)),
      pct_label = sprintf("H %.0f%%\nN %.0f%%", 100 * Heritable, 100 * `Non-heritable`),
      state_label = factor(
        recode(as.character(state), !!!STATE_LABELS),
        levels = rev(unname(STATE_LABELS[STATE_ORDER_FULL]))
      )
    )

  full_limit <- max(abs(rate_diff$log2fc), na.rm = TRUE)
  full_limit <- ceiling(full_limit * 4) / 4

  rate_diff <- rate_diff |>
    mutate(
      fc_label = ifelse(abs(log2fc) >= 0.5, sprintf("%.2f", log2fc), ""),
      text_colour = ifelse(abs(log2fc) >= 1.0, "white", "black")
    )

  abs_rate_limit <- max(rate_abs$frac, na.rm = TRUE)

  fig_B_abs <- ggplot(rate_abs, aes(x = tissue, y = state_label, fill = frac)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    facet_wrap(~ h2_category, nrow = 1) +
    scale_fill_gradientn(
      colours = ABS_RATE_COLORS,
      name = "Overlap",
      limits = c(0, abs_rate_limit),
      breaks = seq(0.1, floor(abs_rate_limit * 10) / 10, by = 0.1),
      labels = percent_format(accuracy = 1)
    ) +
    labs(x = NULL, y = "ChromHMM state") +
    BASE_THEME +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 7.8),
      legend.position = "right",
      legend.key.height = grid::unit(1.0, "cm"),
      legend.key.width = grid::unit(0.35, "cm")
    )

  fig_B_fc <- ggplot(rate_diff, aes(x = tissue, y = state_label, fill = log2fc)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(
      aes(label = fc_label, colour = text_colour),
      size = 2.5,
      show.legend = FALSE
    ) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = EFFECT_COLORS["Higher in non-heritable"],
      mid = "white",
      high = EFFECT_COLORS["Higher in heritable"],
      midpoint = 0,
      limits = c(-full_limit, full_limit),
      oob = squish,
      name = expression(log[2]~"(Her./Non-her.)")
    ) +
    labs(x = NULL, y = "ChromHMM state") +
    BASE_THEME +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 7.8),
      legend.position = "right",
      legend.key.height = grid::unit(1.0, "cm"),
      legend.key.width = grid::unit(0.35, "cm")
    )

  ## Assemble and save

  fig_main <- fig_A / fig_B_quintile +
    plot_layout(heights = c(1.35, 1.0)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 11, face = "bold"))

  fig_hmm_supp <- fig_B_abs / fig_B_fc +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 11, face = "bold"))

  save_plot(fig_main, paste0("repressive_chromatin_main_", pop), width = 10.8, height = 6.2)
  save_plot(fig_A, paste0("repressive_chromatin_fishers_", pop), width = 8.6, height = 5.4)
  save_plot(fig_B_quintile, paste0("repressive_chromatin_quintiles_", pop), width = 7.8, height = 3.5)
  save_plot(fig_hmm_supp, paste0("repressive_chromatin_chromhmm_", pop), width = 7.6, height = 8.6)
  save_plot(fig_A_sens, paste0("repressive_chromatin_sensitivity_peaktype_", pop), width = 7.4, height = 3.8)

}

cat("Figures written to:", OUT_DIR, "\n")

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
if (requireNamespace("sessioninfo", quietly = TRUE)) {
  sessioninfo::session_info()
}
