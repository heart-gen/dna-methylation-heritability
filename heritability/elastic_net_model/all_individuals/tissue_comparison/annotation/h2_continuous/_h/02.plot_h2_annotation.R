#### Visualization: Continuous h² ~ Genomic Annotation ####
##
## Reads results from 01.h2_annotation_glm.R and produces:
##
##   Fig A  — Main: annotation fraction per h² quintile, faceted by region
##            (one colored line per annotation type; 3-column layout)
##   Fig B  — Supplemental: spline-predicted probabilities (continuous h²)
##   Fig C  — Supplemental: odds-ratio forest plot from linear logistic models
##
## Saves PNG + PDF for each figure.
##
## Run from: tissue_comparison/annotation/h2_continuous/_m/
##           conda run -p $ENV Rscript ../_h/02.plot_h2_annotation.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

## Configuration

IN_DIR <- here::here(
  "heritability", "elastic_net_model", "BA_only",
  "tissue_comparison", "annotation", "h2_continuous", "_m"
)
OUT_DIR <- IN_DIR

TISSUE_LABELS <- c(
  "Caudate"      = "Caudate nucleus",
  "DLPFC"        = "DLPFC",
  "Hippocampus"  = "Hippocampus"
)

ANNOT_COLORS <- c(
  "Promoter"           = "#2166AC",   # blue
  "Enhancer"           = "#4DAF4A",   # green
  "1\u20135 kb upstream" = "#984EA3", # purple
  "Intergenic"         = "#D95F02"    # orange
)

ANNOT_ORDER <- c("Promoter", "Enhancer", "1\u20135 kb upstream", "Intergenic")

BASE_THEME <- theme_classic(base_size = 10) +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 10),
    axis.title        = element_text(size = 9),
    axis.text         = element_text(size = 8),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 8),
    legend.key.size   = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    plot.margin       = margin(4, 6, 4, 6)
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

quint  <- fread(file.path(IN_DIR, "quintile_summary.tsv")) |>
  mutate(
    annotation = factor(annotation, levels = ANNOT_ORDER),
    tissue     = factor(tissue,     levels = names(TISSUE_LABELS))
  )

glm    <- fread(file.path(IN_DIR, "glm_linear_results.tsv")) |>
  mutate(
    annotation = factor(annotation, levels = ANNOT_ORDER),
    tissue     = factor(tissue,     levels = names(TISSUE_LABELS)),
    sig        = sig_label(fdr)
  )

spline <- fread(file.path(IN_DIR, "spline_predictions.tsv")) |>
  mutate(
    annotation = factor(annotation, levels = ANNOT_ORDER),
    tissue     = factor(tissue,     levels = names(TISSUE_LABELS))
  )

spline_x_max <- ceiling(max(spline$h2_unscaled, na.rm = TRUE) * 10) / 10

## Figure A — Quintile fraction plot (main figure)

fig_a <- ggplot(quint,
    aes(x = h2_quintile, y = frac,
        color = annotation, group = annotation)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = annotation),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8, shape = 19) +
  facet_wrap(
    ~ tissue,
    nrow   = 1,
    labeller = as_labeller(TISSUE_LABELS)
  ) +
  scale_color_manual(values = ANNOT_COLORS, breaks = ANNOT_ORDER) +
  scale_fill_manual( values = ANNOT_COLORS, breaks = ANNOT_ORDER) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = expression(paste("h"^2, " quintile")),
    y = "VMRs overlapping annotation (%)"
  ) +
  BASE_THEME +
  theme(
    legend.position  = "right",
    legend.direction = "vertical"
  )

save_plot(fig_a, "h2_annotation_quintile", width = 7.2, height = 3.2)

## Figure B — Spline-predicted probabilities (supplemental)

fig_b <- ggplot(spline,
    aes(x = h2_unscaled, y = prob,
        color = annotation, fill = annotation)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.13, color = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(
    ~ tissue,
    nrow     = 1,
    labeller = as_labeller(TISSUE_LABELS)
  ) +
  scale_color_manual(values = ANNOT_COLORS, breaks = ANNOT_ORDER) +
  scale_fill_manual( values = ANNOT_COLORS, breaks = ANNOT_ORDER) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 5),
    limits = c(0, spline_x_max),
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = expression(paste("SNP-based heritability (h"^2, ")")),
    y = "Predicted probability of\nannotation overlap (%)"
  ) +
  BASE_THEME +
  theme(legend.position = "right")

save_plot(fig_b, "h2_annotation_spline", width = 7.2, height = 3.2)

## Figure C — OR forest plot (supplemental)

glm_plot <- glm |>
  mutate(
    tissue_label = recode(as.character(tissue), !!!TISSUE_LABELS),
    label_right  = sprintf("%.2f [%.2f\u2013%.2f]%s",
                            estimate, conf.low, conf.high, sig)
  )

# Nudge significance stars slightly above error bars
fig_c <- ggplot(glm_plot,
    aes(x = estimate, y = annotation,
        color = annotation, xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60",
             linewidth = 0.5) +
  geom_errorbarh(height = 0.2, linewidth = 0.6) +
  geom_point(size = 2.4) +
  geom_text(aes(x = conf.high, label = sig),
            hjust = -0.3, vjust = 0.5, size = 3, color = "black") +
  facet_wrap(
    ~ tissue_label,
    nrow     = 1,
    scales   = "free_x"
  ) +
  scale_color_manual(values = ANNOT_COLORS, breaks = ANNOT_ORDER,
                     guide = "none") +
  scale_x_continuous(
    trans  = "log2",
    labels = number_format(accuracy = 0.01)
  ) +
  labs(
    x = expression(paste(
      "Odds ratio per 0.1 unit increase in h"^2, " (log"[2], " scale)"
    )),
    y = NULL
  ) +
  BASE_THEME +
  theme(
    axis.text.y  = element_text(size = 8),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.major.y = element_blank()
  )

save_plot(fig_c, "h2_annotation_OR", width = 7.2, height = 3.0)

## Combined supplemental figure (B + C stacked)

fig_supp <- fig_b / fig_c +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))

save_plot(fig_supp, "h2_annotation_supplemental", width = 7.2, height = 6.4)

## Sensitivity: all VMRs (including low-prediction)

quint_sa <- fread(file.path(IN_DIR, "quintile_summary_sensitivity.tsv")) |>
  mutate(
    annotation = factor(annotation, levels = ANNOT_ORDER),
    tissue     = factor(tissue,     levels = names(TISSUE_LABELS))
  )

fig_sa <- ggplot(quint_sa,
    aes(x = h2_quintile, y = frac,
        color = annotation, group = annotation)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = annotation),
              alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_wrap(~ tissue, nrow = 1,
             labeller = as_labeller(TISSUE_LABELS)) +
  scale_color_manual(values = ANNOT_COLORS, breaks = ANNOT_ORDER) +
  scale_fill_manual( values = ANNOT_COLORS, breaks = ANNOT_ORDER) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(
    x = expression(paste("h"^2, " quintile (all VMRs)")),
    y = "VMRs overlapping annotation (%)"
  ) +
  BASE_THEME +
  theme(legend.position = "right")

save_plot(fig_sa, "h2_annotation_sensitivity", width = 7.2, height = 3.2)

cat("Figures written to:", OUT_DIR, "\n")

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
