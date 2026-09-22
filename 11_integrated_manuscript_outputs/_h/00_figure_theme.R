#### 11_integrated_manuscript_outputs / 00_figure_theme: shared figure grammar ####
##
## One source for the manuscript's visual language. The v1 tree carried this
## same theme, palette, and save helper copy-pasted into ~40 scripts across
## three cohort trees, which is how panels drifted apart between figures. The
## definitions here are ported from
##   local-snp-prediction/all_individuals/tissue_comparison/annotation/
##     repressive_chromatin/_h/02.plot_repressive_chromatin.R  (lines 88-131)
## so new figures sit alongside the existing ones without a style break.
##
## Not a standalone script: sourced by the module 11 figure and table builders.
##
## WHY THIS LIVES IN _h/ AND NOT 00_shared/
## `_m/runs/*/code/` snapshots `_h/` and `config/`, but NOT `00_shared/`. A
## defect in shared code therefore reaches a sealed run with no provenance
## trail -- which is exactly how a `block_jackknife_cov()` bug reached three
## sealed Module 10 runs on 2026-09-20. Figure style stays inside the snapshot.
##
## v1 STYLE DECISIONS THIS FILE SETTLES
## v1 was internally inconsistent about brain-region colour: two scripts assign
## Caudate #B36F61 / DLPFC #7372A6 / Hippocampus #E3C962, four assign the first
## two swapped. Rather than revive an inconsistency, v2 keeps the Era B
## repressive-chromatin register below and states the choice once, here.
## v1's Heritable / Non-heritable / Low-prediction palette is deliberately NOT
## ported: those classes are retired (AGENTS.md 3).

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})

## ------------------------------------------------------------------- theme
##
## base_size 10 with axis text at 8: sized for a full-width (180 mm) panel
## reduced into a two-column PDF, where 8 pt is about the floor for legibility.
##
## The font family is named rather than left to the device default, so a PDF
## rendered on a compute node and one rendered on a laptop embed the same
## faces. "sans" resolves to Helvetica in the pdf device and to the system sans
## face under cairo; change it in this one place if the journal asks.

FIG_FONT <- "sans"

BASE_THEME <- theme_classic(base_size = 10, base_family = FIG_FONT) +
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

## v1's Era B value-axis gridline, which several v1 scripts add and v2 was
## pasting inline. Add it to a panel whose reader has to compare magnitudes
## across a category axis; leave it off dense line and density panels.
GRID_Y <- theme(panel.grid.major.y = element_line(colour = "grey92",
                                                  linewidth = 0.3))
GRID_X <- theme(panel.grid.major.x = element_line(colour = "grey92",
                                                  linewidth = 0.3))

## ------------------------------------------------------------- panel tags
##
## LOWERCASE. `content/91.figure-legends.md` cites panels as **a.**, **b.**,
## and v1's hand-assembled SVGs were relabelled a/b/c by hand after being
## rendered with tag_levels = "A". Generating them lowercase is what removes
## that manual step. Bold, 11 pt, matching v1
## (02.plot_h2_annotation.R:220-222).

TAG_THEME <- theme(plot.tag = element_text(face = "bold", size = 11))

#' patchwork annotation for a multi-panel figure.
#'
#' Use this instead of hand-writing plot_annotation(), so no builder can
#' reintroduce uppercase tags. Combine with `& TAG_THEME`:
#'
#'     figure <- (pA / pB) + fig_tags() & TAG_THEME
fig_tags <- function(...) {
    patchwork::plot_annotation(tag_levels = "a", ...)
}

## ----------------------------------------------------------------- palette
##
## The v1 warm register. Ordered so the first three stay distinguishable in
## grayscale and under deuteranopia; grey is reserved for null/not-significant
## AND for an estimate deliberately set aside from a claim (see Module 04's
## technically-confounded caudate LINE/L1).

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

## --------------------------------------------------------- upstream run IDs
##
## Module 01's QC-refresh runs carry array coverage and genomic context, which
## the accepted catalog runs predate. They are a refresh of an accepted run
## rather than an accepted run of their own, so `require_accepted_upstream()`
## cannot resolve them and they must be named. Named ONCE, here, because two
## builders read them and a stale copy in one of them is how Figure 2 came to
## cite a retired Module 02 run.
##
## Every other upstream run ID is resolved through
## `00_shared/gates.R::require_accepted_upstream()` against the README of
## record. Do not add run-ID string templates to a builder.

QC_REFRESH_RUN <- function(cohort, region) {
    paste0("vmrcatqc-", cohort, "-", region, "-20260826-a")
}

## -------------------------------------------------------------- statistics

#' v1's significance key, identical across every v1 plotting script.
#'
#' Legends spell it out as P < 0.05 (*), P < 0.01 (**), P < 0.001 (***).
#' Feed this FDR-adjusted values; the manuscript reports q, not nominal p.
sig_stars <- function(q) {
    data.table::fifelse(
        is.na(q), "",
        data.table::fifelse(q < 0.001, "***",
        data.table::fifelse(q < 0.01, "**",
        data.table::fifelse(q < 0.05, "*", ""))))
}

SIG_KEY <- "*** q < 0.001   ** q < 0.01   * q < 0.05"

#' v1's diverging log2 odds-ratio fill, ported from
#' annotation/enrichment/_h/02.plot_heatmap.R:76-84.
#'
#' The symmetric limit is taken from the 90th percentile of |log2 OR| with a
#' floor of 0.3, and out-of-bounds values are squished rather than dropped, so
#' a single extreme cell cannot wash the panel out. Pass the values that will
#' be plotted.
scale_fill_log2or <- function(values, name = expression(log[2]~"(OR)"), ...) {
    v <- abs(values[is.finite(values)])
    lim <- max(ceiling(stats::quantile(v, 0.90, na.rm = TRUE) * 10) / 10, 0.3)
    ggplot2::scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
        limits = c(-lim, lim), oob = scales::squish, name = name,
        breaks = scales::pretty_breaks(n = 5), ...)
}

#' Horizontal error bar.
#'
#' geom_errorbarh() is deprecated as of ggplot2 4.0.0 (this environment runs
#' 4.0.1) in favour of geom_errorbar(orientation = "y"). Wrapped so the
#' deprecation is handled in one place and the `height` -> `width` rename does
#' not have to be remembered at every call site.
errorbar_h <- function(..., width = 0) {
    ggplot2::geom_errorbar(..., orientation = "y", width = width)
}

## -------------------------------------------------------------------- output

FIG_WIDTH_FULL   <- 7.09   # 180 mm
FIG_WIDTH_THREEQ <- 5.51   # 140 mm
FIG_WIDTH_SINGLE <- 3.46   # 88 mm

## A main figure has to print on one journal page. v1's Era B figures ran
## 2.5-6.4 in tall; the first v2 drafts reached 11.4 in, which no layout
## accommodates. save_figure() refuses anything taller rather than leaving it
## to a visual check nobody repeats.
FIG_HEIGHT_MAX <- 9.5

## PNG is a review copy, not the deliverable: 300 dpi matches v1 Era B.
FIG_DPI <- 300

#' Write a figure as PDF (vector, for the journal), PNG (review copies) and
#' SVG (drop-in for the manubot manuscript build, which embeds .svg/.png).
#'
#' THE PDF DEVICE IS cairo_pdf, NOT pdf. R's base pdf() device writes
#' single-byte encodings and silently drops any character it cannot map. The
#' axis label "Out-of-fold \u03c1\u00b2" in Figure 2 panel a triggered exactly
#' that:
#'   conversion failure on 'Out-of-fold rho-squared' in 'mbcsToSbcs'
#' -- a WARNING, in a SLURM log, on the file that goes to the journal, while
#' the PNG review copy rendered correctly. Greek letters, superscripts, en
#' dashes and multiplication signs are unavoidable in these axis labels, so the
#' device has to be UTF-8 capable. cairo_pdf is, and embeds its fonts.
#'
#' bg = "white" is set on every device. Without it ggsave writes a transparent
#' background, which reads as a black panel in several PDF viewers and in
#' anything that composites on a dark surface.
#'
#' svglite is not in the epigenomics environment, so SVG goes through
#' grDevices::svg(), also cairo. It renders glyphs as paths, so SVG text is not
#' selectable -- fine for embedding, which is all the manuscript build needs.
#'
#' If cairo is unavailable, both fall back to the base devices with an explicit
#' warning rather than silently shipping mangled labels.
save_figure <- function(plot_obj, name, width, height, fig_dir, svg = TRUE) {
    if (height > FIG_HEIGHT_MAX) {
        stop("Figure '", name, "' is ", height, " in tall; the limit is ",
             FIG_HEIGHT_MAX, " in (one journal page). Move panels to a ",
             "supplementary figure rather than raising the limit.")
    }
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    has_cairo <- isTRUE(capabilities("cairo"))
    if (!has_cairo) {
        warning("No cairo support: PDF falls back to pdf(), which drops ",
                "non-Latin-1 glyphs such as Greek letters. Check every axis ",
                "label in ", name, " before using the PDF.")
    }

    ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot = plot_obj,
           width = width, height = height, bg = "white",
           device = if (has_cairo) grDevices::cairo_pdf else grDevices::pdf)
    ggsave(file.path(fig_dir, paste0(name, ".png")), plot = plot_obj,
           width = width, height = height, dpi = FIG_DPI, bg = "white")

    if (svg && has_cairo) {
        f <- file.path(fig_dir, paste0(name, ".svg"))
        grDevices::svg(filename = f, width = width, height = height,
                       bg = "white")
        on.exit(grDevices::dev.off(), add = TRUE)
        print(plot_obj)
    } else if (svg) {
        warning("No cairo support; skipping SVG for ", name,
                " (PDF and PNG were written).")
    }

    message("[figure] ", name, " (", width, " x ", height, " in)")
    invisible(file.path(fig_dir, name))
}

## DENSE LAYERS
## v1 shipped supplementary SVGs of 9.4 MB and 14.7 MB because per-VMR scatter
## points were written as vector geometry. ggrastr is not installed and is not
## worth a new dependency: summarise dense per-VMR layers to deciles, bins or
## hexes before plotting, as Figures 1 and 2 already do. If a genuinely dense
## layer is unavoidable, say so in the panel's source data.

#' Write a panel's source data with the provenance AGENTS.md 7.11 requires:
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

#' Name a panel's source data by the tag patchwork will RENDER, not by the R
#' variable that holds the plot.
#'
#' Figure 1 assembled (pA / pB / (pC+pD) / (pF+pG) / pE), so pF rendered as
#' "e", pG as "f" and pE as "g" -- while its sd() calls named them panelF,
#' panelG and panelE. The source data disagreed with the figure it documented.
#' Pass the plots in assembly order and read the tag off the position.
panel_tags <- function(n) letters[seq_len(n)]
