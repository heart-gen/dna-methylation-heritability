#### Visualization: Intergenic VMR Open Chromatin & ABC Enrichment ####
##
## Produces manuscript-quality figures from scripts 01 and 02:
##
##   Fig 1  — Cell-type enrichment heatmap (Fisher's OR, heritable vs non-heritable)
##   Fig 2  — Forest plot: OR + 95% CI per cell type (per tissue + pooled)
##   Fig 3  — Quintile plot: open chromatin overlap fraction vs h² quintile
##   Fig 4  — Cell-type specificity: stacked bar of n_celltypes_overlap
##   Fig 5  — GO pathway dot plot (heritable vs non-heritable target genes)
##   Fig 6  — ABC score distribution by h2 category (violin + boxplot)
##
## Main figure: Fig 1 + Fig 3 (left-to-right narrative: enrichment → gradient)
## Supplemental: Fig 2 + Fig 4 + Fig 5 + Fig 6
##
## Run: conda run -p $ENV_PATH/epigenomics Rscript ../_h/03.plot_results.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(RColorBrewer)
})

## Configuration

IN_DIR <- here::here(
  "heritability", "elastic_net_model", "all_individuals",
  "tissue_comparison", "annotation", "open_chromatin", "_m"
)
OUT_DIR <- IN_DIR

POPULATIONS <- c("AA", "EA")

TISSUE_LABELS <- c(
  "Caudate"     = "Caudate nucleus",
  "DLPFC"       = "DLPFC",
  "Hippocampus" = "Hippocampus",
  "Pooled"      = "Pooled"
)
TISSUE_ORDER <- c("Caudate", "DLPFC", "Hippocampus", "Pooled")

# Cell-type display labels and order (neuronal first, then glial, then other)
CT_LABELS <- c(
  "Union" = "All cell types",
  "Exc"   = "Excitatory",
  "Inh"   = "Inhibitory",
  "Astro" = "Astrocyte",
  "Oligo" = "Oligodendrocyte",
  "OPC"   = "OPC",
  "Micro" = "Microglia",
  "Endo"  = "Endothelial"
)
CT_ORDER <- names(CT_LABELS)  # Union last in heatmap; first in forest plot

# Accessible color palette
H2_COLORS <- c(
  "Heritable"     = "#1B7837",
  "Non-heritable" = "#762A83"
)
# Colors keyed by display label (must match CT_LABELS values)
CT_COLORS <- c(
  "All cell types"  = "#555555",
  "Excitatory"      = "#E41A1C",
  "Inhibitory"      = "#377EB8",
  "Astrocyte"       = "#FF7F00",
  "Oligodendrocyte" = "#4DAF4A",
  "OPC"             = "#984EA3",
  "Microglia"       = "#A65628",
  "Endothelial"     = "#F781BF"
)

BASE_THEME <- theme_classic(base_size = 10) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 10),
    axis.title         = element_text(size = 9),
    axis.text          = element_text(size = 8),
    legend.title       = element_text(size = 8),
    legend.text        = element_text(size = 8),
    legend.key.size    = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    plot.margin        = margin(4, 6, 4, 6)
  )

save_plot <- function(p, name, width, height) {
  for (ext in c(".pdf", ".png")) {
    ggsave(
      filename = file.path(OUT_DIR, paste0(name, ext)),
      plot     = p,
      width    = width,
      height   = height,
      dpi      = 300
    )
  }
  invisible(p)
}

sig_label <- function(fdr) {
  dplyr::case_when(
    fdr < 0.001 ~ "***",
    fdr < 0.01  ~ "**",
    fdr < 0.05  ~ "*",
    TRUE        ~ ""
  )
}

## Load data

for (pop in POPULATIONS) {

  cat(sprintf("Processing population: %s \n", pop))

  fisher_df  <- fread(file.path(IN_DIR, paste0("fishers_heritable_vs_nonheritable_", pop, ".tsv")))
  logistic_df <- fread(file.path(IN_DIR, paste0("logistic_continuous_h2_", pop, ".tsv")))
  quintile_df <- fread(file.path(IN_DIR, paste0("quintile_open_chromatin_", pop, ".tsv")))
  vmr_df      <- fread(file.path(IN_DIR, paste0("intergenic_vmr_atac_overlap_", pop, ".tsv")))

  go_file <- file.path(IN_DIR, paste0("enrichr_go_bp_all_", pop, ".tsv"))
  go_df   <- if (file.exists(go_file)) fread(go_file) else NULL

  abc_file <- file.path(IN_DIR, paste0("abc_vmr_gene_links_", pop, ".tsv"))
  abc_df   <- if (file.exists(abc_file)) fread(abc_file) else NULL

  # Factor ordering
  fisher_df <- fisher_df |>
    mutate(
      cell_type = factor(CT_LABELS[cell_type], levels = CT_LABELS),
      tissue    = factor(tissue, levels = TISSUE_ORDER),
      sig       = sig_label(fdr)
    )

  quintile_df <- quintile_df |>
    mutate(
      cell_type = factor(CT_LABELS[cell_type], levels = CT_LABELS),
      tissue    = factor(tissue, levels = TISSUE_ORDER)
    )

  logistic_df <- logistic_df |>
    mutate(
      cell_type = factor(CT_LABELS[cell_type], levels = CT_LABELS),
      tissue    = factor(tissue, levels = TISSUE_ORDER),
      sig       = sig_label(fdr)
    )

  ## Fig 1 — Cell-type enrichment tile heatmap
  ## OR (heritable vs non-heritable intergenic VMRs) × cell type × tissue

  # Clip extreme ORs for display
  or_lim <- c(0.5, 3.0)
  heatmap_df <- fisher_df |>
    mutate(
      or_clipped = pmin(pmax(or, or_lim[1]), or_lim[2]),
      log2_or    = log2(or_clipped)
    )

  lim_log2 <- log2(or_lim)

  fig1 <- ggplot(heatmap_df,
      aes(x = tissue, y = cell_type, fill = log2_or)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sig), size = 3, vjust = 0.75, color = "white") +
    scale_fill_gradient2(
      low      = "#762A83",
      mid      = "white",
      high     = "#1B7837",
      midpoint = 0,
      limits   = lim_log2,
      oob      = scales::squish,
      name     = expression(log[2]~OR~"(Her. / Non-her.)"),
      labels   = function(x) sprintf("%.1f", x)
    ) +
    scale_x_discrete(labels = TISSUE_LABELS) +
    labs(x = NULL, y = NULL) +
    BASE_THEME +
    theme(
      axis.text.x        = element_text(angle = 30, hjust = 1, size = 9),
      axis.text.y        = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      legend.position    = "right",
      legend.key.height  = unit(1.2, "cm"),
      legend.key.width   = unit(0.35, "cm")
    )

  save_plot(fig1, paste0("fig1_celltype_heatmap_", pop), width = 5.0, height = 3.5)

  ## Fig 2 — Forest plot: OR + 95% CI per cell type, faceted by tissue (supplemental)

  fig2 <- ggplot(fisher_df,
      aes(x = or, y = cell_type, color = cell_type,
          xmin = ci_lo, xmax = ci_hi)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60",
              linewidth = 0.4) +
    geom_errorbarh(height = 0.25, linewidth = 0.6) +
    geom_point(size = 2.2) +
    geom_text(aes(x = ci_hi, label = sig),
              hjust = -0.3, vjust = 0.5, size = 3, color = "black") +
    facet_wrap(~ tissue, nrow = 1, labeller = as_labeller(TISSUE_LABELS),
              scales = "free_x") +
    scale_color_manual(values = CT_COLORS, guide = "none") +
    scale_x_continuous(trans  = "log2",
                      labels = number_format(accuracy = 0.01)) +
    labs(x = expression("Odds ratio (heritable / non-heritable), log"[2]~scale),
        y = NULL) +
    BASE_THEME +
    theme(
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      axis.text.y        = element_text(size = 8)
    )

  save_plot(fig2, paste0("fig2_forest_plot_", pop), width = 8.0, height = 3.5)

  ## Fig 3 — Quintile plot: fraction overlapping open chromatin vs h² quintile

  # Show selected cell types for main figure clarity
  main_cts <- CT_LABELS[c("Union", "Exc", "Inh", "Astro", "Oligo", "Micro")]

  quint_main <- quintile_df |>
    filter(
      tissue %in% c("DLPFC", "Caudate", "Hippocampus"),
      cell_type %in% main_cts
    ) |>
    mutate(cell_type = droplevels(cell_type))

  ct_colors_main <- CT_COLORS[levels(quint_main$cell_type)]

  fig3 <- ggplot(quint_main,
      aes(x = h2_quintile, y = frac,
          color = cell_type, group = cell_type)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = cell_type),
                alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8, shape = 19) +
    facet_wrap(~ tissue, nrow = 1, labeller = as_labeller(TISSUE_LABELS)) +
    scale_color_manual(values = ct_colors_main, name = "Cell type") +
    scale_fill_manual( values = ct_colors_main, name = "Cell type") +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.06))
    ) +
    labs(
      x = expression(paste("h"^2, " quintile (intergenic VMRs)")),
      y = "VMRs overlapping open chromatin (%)"
    ) +
    BASE_THEME +
    theme(legend.position = "right")

  save_plot(fig3, paste0("fig3_quintile_open_chromatin_", pop), width = 7.5, height = 3.2)

  ## Fig 4 — Cell-type specificity stacked bar (supplemental)
  ## How many cell types does each heritable intergenic VMR overlap?

  ct_spec_df <- vmr_df |>
    filter(h2_category %in% c("Heritable", "Non-heritable")) |>
    mutate(
      n_cts_label = factor(
        pmin(n_celltypes_overlap, 5),
        levels = 0:5,
        labels = c("0", "1", "2", "3", "4", "5+")
      ),
      h2_category = factor(h2_category, levels = c("Heritable", "Non-heritable"))
    ) |>
    count(tissue, h2_category, n_cts_label) |>
    group_by(tissue, h2_category) |>
    mutate(frac = n / sum(n)) |>
    ungroup()

  fig4 <- ggplot(ct_spec_df,
      aes(x = h2_category, y = frac, fill = n_cts_label)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    facet_wrap(~ tissue, nrow = 1, labeller = as_labeller(TISSUE_LABELS)) +
    scale_fill_manual(
      values = RColorBrewer::brewer.pal(6, "YlOrRd"),
      name   = "Cell types\noverlapping"
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                      expand = expansion(mult = c(0, 0.03))) +
    scale_x_discrete(labels = c("Heritable" = "Her.", "Non-heritable" = "Non-her.")) +
    labs(x = NULL, y = "Fraction of intergenic VMRs") +
    BASE_THEME +
    theme(
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      legend.position    = "right"
    )

  save_plot(fig4, paste0("fig4_celltype_specificity_", pop), width = 6.5, height = 3.2)

  ## Fig 5 — GO pathway dot plot (supplemental)

  if (!is.null(go_df) && nrow(go_df) > 0) {

    COMPARISON_LABELS <- c(
      "Pooled_Heritable"    = "Heritable",
      "Pooled_NonHeritable" = "Non-heritable",
      "Pooled_Q5"           = "Top h² (Q5)",
      "Pooled_Q1Q4"         = "Lower h² (Q1–Q4)"
    )

    go_plot_df <- go_df |>
      filter(comparison %in% names(COMPARISON_LABELS)) |>
      mutate(
        comparison = factor(COMPARISON_LABELS[comparison],
                            levels = COMPARISON_LABELS),
        log10_fdr  = -log10(p.adjust),
        GeneRatio_num = sapply(GeneRatio, function(x) {
          parts <- strsplit(x, "/")[[1]]
          as.numeric(parts[1]) / as.numeric(parts[2])
        })
      ) |>
      arrange(p.adjust) |>
      # Top 10 per comparison
      group_by(comparison) |>
      slice_min(p.adjust, n = 10) |>
      ungroup()

    # Keep only terms appearing in at least one comparison (de-duplicate Description)
    top_terms <- go_plot_df |>
      distinct(comparison, Description) |>
      count(Description) |>
      pull(Description)

    go_plot_df <- go_plot_df |>
      filter(Description %in% top_terms) |>
      mutate(Description = factor(Description,
        levels = rev(unique(go_plot_df$Description[order(go_plot_df$p.adjust)]))))

    fig5 <- ggplot(go_plot_df,
        aes(x = comparison, y = Description,
            size = GeneRatio_num, color = log10_fdr)) +
      geom_point() +
      scale_color_gradient(
        low    = "#FEE08B",
        high   = "#D73027",
        name   = expression(-log[10]~FDR)
      ) +
      scale_size_continuous(
        range = c(1.5, 5),
        name  = "Gene ratio"
      ) +
      labs(x = NULL, y = NULL) +
      BASE_THEME +
      theme(
        axis.text.x  = element_text(angle = 30, hjust = 1, size = 8),
        axis.text.y  = element_text(size = 7.5),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
      )

    save_plot(fig5, paste0("fig5_go_dotplot_", pop), width = 6.5, height = 5.5)
  } else {
    message("GO enrichment results not found — skipping Fig 5")
  }

  ## Fig 6 — ABC score distribution by h2 category (supplemental)

  if (!is.null(abc_df) && nrow(abc_df) > 0) {
    abc_viol_df <- abc_df |>
      filter(tissue == "Pooled",
            h2_category %in% c("Heritable", "Non-heritable")) |>
      mutate(h2_category = factor(h2_category,
                                  levels = c("Heritable", "Non-heritable")))

    wt <- wilcox.test(ABC.Score ~ h2_category, data = abc_viol_df)
    wt_label <- sprintf("Wilcoxon p = %.2g", wt$p.value)

    y_max <- quantile(abc_viol_df$ABC.Score, 0.99, na.rm = TRUE)

    fig6 <- ggplot(abc_viol_df,
        aes(x = h2_category, y = ABC.Score, fill = h2_category)) +
      geom_violin(trim = TRUE, alpha = 0.75, color = NA, scale = "width") +
      geom_boxplot(width = 0.12, fill = "white", color = "grey30",
                  outlier.shape = NA, linewidth = 0.5) +
      annotate("text", x = 1.5, y = y_max * 1.05,
              label = wt_label, size = 3, color = "grey20") +
      scale_fill_manual(values = H2_COLORS, guide = "none") +
      scale_y_continuous(
        limits = c(0, y_max * 1.1),
        expand = expansion(mult = c(0, 0.02))
      ) +
      scale_x_discrete(labels = c("Heritable" = "Heritable\nintergenic",
                                  "Non-heritable" = "Non-heritable\nintergenic")) +
      labs(x = NULL, y = "ABC score") +
      BASE_THEME +
      theme(panel.grid.major.x = element_blank())

    save_plot(fig6, paste0("fig6_abc_score_violin_", pop), width = 3.2, height = 3.5)
  } else {
    message("ABC links file not found — skipping Fig 6")
  }

  ## Main figure: Fig 1 (heatmap) | Fig 3 (quintile) — side by side

  fig_main <- fig1 + fig3 +
    plot_layout(widths = c(1, 1.8)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))

  save_plot(fig_main, paste0("main_fig_cre_enrichment_", pop), width = 12.5, height = 3.8)

  ## Supplemental figure: Fig 2 + Fig 4 stacked

  supp_top <- fig2 + fig4 +
    plot_layout(ncol = 2, widths = c(1.3, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))

  save_plot(supp_top, paste0("supp_fig_forest_celltype_", pop), width = 14.5, height = 3.8)

}

cat("All figures saved to:", OUT_DIR, "\n")

#### Reproducibility ####
cat("\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
