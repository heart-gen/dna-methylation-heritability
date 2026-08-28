#### 10_integrated_manuscript_outputs / 00_figure_theme: shared figure grammar ####
##
## One source for the manuscript's visual language. The v1 tree carried this
## same theme, palette, and save helper copy-pasted into ~40 scripts across
## three cohort trees, which is how panels drifted apart between figures. The
## definitions here are ported from
##   local-snp-prediction/all_individuals/tissue_comparison/annotation/
##     repressive_chromatin/_h/02.plot_repressive_chromatin.R  (lines 88-131)
## so new figures sit alongside the existing ones without a style break.
##
## Not a standalone script: sourced by the module 10 figure builders.

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})

## ------------------------------------------------------------------- theme
##
## base_size 10 with axis text at 8: sized for a full-width (180 mm) panel
## reduced into a two-column PDF, where 8 pt is about the floor for legibility.

BASE_THEME <- theme_classic(base_size = 10) +
    theme(
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 10),
        axis.title       = element_text(size = 9),
        axis.text        = element_text(size = 8, colour = "black"),
        legend.title     = element_text(size = 8),
        legend.text      = element_text(size = 8),
        legend.key.size  = grid::unit(0.42, "cm"),
        panel.border     = element_blank(),
        plot.margin      = margin(5, 8, 5, 8)
    )

## Figure-level interpretation belongs in the caption, never inside the panel
## (AGENTS.md 11 and the manuscript-figures convention).
NO_TITLES <- theme(plot.title = element_blank(), plot.subtitle = element_blank())

## ----------------------------------------------------------------- palette
##
## The v1 warm register. Ordered so the first three stay distinguishable in
## grayscale and under deuteranopia; grey is reserved for null/not-significant.

PAL_RUST     <- "#B6523A"
PAL_BLUE     <- "#3A6F8F"
PAL_BLUE_LT  <- "#5B8FA8"
PAL_GREEN    <- "#6FA287"
PAL_TAN      <- "#C17C59"
PAL_CHARCOAL <- "#3F3F3F"
PAL_NULL     <- "#B6B6B6"

REGION_LABELS <- c(caudate = "Caudate", dlpfc = "DLPFC",
                   hippocampus = "Hippocampus")
REGION_ORDER  <- unname(REGION_LABELS)
REGION_COLORS <- c(Caudate = PAL_RUST, DLPFC = PAL_BLUE,
                   Hippocampus = PAL_GREEN)

COHORT_LABELS <- c(AA = "Black American donors",
                   all_individuals = "All donors")

## Sequential ramp for heatmap fills, light to dark.
SEQ_RAMP <- c("#F7F3EE", "#D49A72", "#8E4426")

PLATFORM_COLORS <- c("450K" = PAL_TAN, "EPIC" = PAL_CHARCOAL)

as_region <- function(x) factor(unname(REGION_LABELS[x]), levels = REGION_ORDER)

## -------------------------------------------------------------------- output

FIG_WIDTH_FULL   <- 7.09   # 180 mm
FIG_WIDTH_THREEQ <- 5.51   # 140 mm
FIG_WIDTH_SINGLE <- 3.46   # 88 mm

#' Write a figure as PDF (vector, for the journal) and PNG (for review copies).
save_figure <- function(plot_obj, name, width, height, fig_dir) {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot = plot_obj,
           width = width, height = height, useDingbats = FALSE)
    ggsave(file.path(fig_dir, paste0(name, ".png")), plot = plot_obj,
           width = width, height = height, dpi = 400)
    message("[figure] ", name, " (", width, " x ", height, " in)")
    invisible(file.path(fig_dir, name))
}

#' Write a panel's source data with the provenance AGENTS.md 7.9 requires:
#' every panel must record its source run ID, table, script, and filter.
write_source_data <- function(dt, name, source_run_id, source_table, script,
                              filter_desc, data_dir) {
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    out <- data.table::as.data.table(dt)
    out[, `:=`(source_run_id = paste(unique(source_run_id), collapse = ";"),
               source_table  = source_table,
               source_script = script,
               row_filter    = filter_desc)]
    write_atomic(out, file.path(data_dir, paste0(name, ".tsv")))
    invisible(out)
}
