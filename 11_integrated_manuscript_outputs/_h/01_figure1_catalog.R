#### 10 / Figure 1: cohort, corrected VMR catalog, and off-array coverage ####
##
## AGENTS.md 11 assigns Figure 1 the cohort, the corrected VMR catalog, and WGBS
## coverage outside array-accessible CpGs. AGENTS.md 7.9 makes this module the
## only place manuscript figures may be assembled.
##
## Panels, in the order patchwork renders them -- the tag a figure SHOWS is
## the only panel name this script may use, including in source_data/.
##   a  donors and CpGs assayed per region
##   b  VMRs per chromosome, three regions
##   c  VMR width and CpGs per VMR
##   d  WGBS coverage outside array-accessible CpGs   <- the 2.2 contribution
##   e  genomic compartment the VMRs occupy (descriptive, not enrichment)
##   f  distance to the nearest gene
##
## Catalog turnover against the invalid legacy calls is a repair audit against
## a catalog AGENTS.md 8 calls invalid, so it ships as its own supplementary
## figure rather than as a seventh panel here. Dropping it is also what brings
## the main figure under one journal page.
##
## Usage:
##   Rscript 01_figure1_catalog.R --cohort AA --run-id fig-all-20260826
##   Rscript 01_figure1_catalog.R --cohort AA --run-id ID --platform EPIC

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(V2_ROOT, "11_integrated_manuscript_outputs", "_h",
                 "00_figure_theme.R"))

suppressPackageStartupMessages({
    library(patchwork)
    library(scales)
})

opts <- parse_v2_args(require = c("cohort", "run_id"))
cohort <- opts$cohort
## 450K in the main figure; EPIC is the stricter supplemental comparator.
platform <- if (is.null(opts$platform)) "450K" else toupper(opts$platform)
regions <- load_config("cohorts")$regions

module_root <- file.path(V2_ROOT, "11_integrated_manuscript_outputs")
## --out-dir renders a review draft outside the immutable run tree, so a panel
## can be iterated on without minting and sealing a run. Same flag as
## 06_figure_region_donor_generalization.R.
run_dir  <- if (!is.null(opts$out_dir)) opts$out_dir else
    file.path(module_root, "_m", "runs", opts$run_id)
fig_dir  <- file.path(run_dir, "figures")
data_dir <- file.path(run_dir, "source_data")

## The catalog run holds the VMR calls; the qc refresh holds array coverage and
## genomic context, which the sealed catalog runs predate.
##
## The catalog run is resolved through the acceptance gate (AGENTS.md 6), not
## from a run-ID string template. A template is how Figure 2 came to cite a
## Module 02 run that had been retired months earlier: nothing fails when the
## README of record moves on, the figure just silently keeps building on a
## superseded run. The QC refresh has no acceptance row of its own -- it
## re-runs QC over an already accepted catalog -- so it is named, once, in
## 00_figure_theme.R.
CATALOG_RUN <- local({
    cache <- new.env(parent = emptyenv())
    function(r) {
        if (is.null(cache[[r]])) {
            cache[[r]] <- require_accepted_upstream("01_vmr_catalog",
                                                    cohort = cohort,
                                                    region = r)$run_id
        }
        cache[[r]]
    }
})
QC_RUN <- function(r) QC_REFRESH_RUN(cohort, r)

cat_dir <- function(r) file.path(V2_ROOT, "01_vmr_catalog", "_m", "runs", CATALOG_RUN(r))
qc_dir  <- function(r) file.path(V2_ROOT, "01_vmr_catalog", "_m", "runs", QC_RUN(r))

read_by_region <- function(fun) {
    rbindlist(lapply(regions, function(r) {
        d <- fun(r); if (is.null(d) || nrow(d) == 0) return(NULL)
        d[, region := r][]
    }), fill = TRUE)
}

vmr      <- read_by_region(function(r) fread(file.path(cat_dir(r), "vmr", "vmr_catalog.tsv")))
cutoffs  <- read_by_region(function(r) fread(file.path(cat_dir(r), "vmr", "sd_cutoffs.tsv")))
tqc      <- read_by_region(function(r) fread(file.path(qc_dir(r), "qc", "technical_qc.tsv")))
coverage <- read_by_region(function(r) fread(file.path(qc_dir(r), "qc", "array_coverage.tsv")))
turnover <- read_by_region(function(r) fread(file.path(qc_dir(r), "qc", "vmr_turnover.tsv")))
gctx     <- read_by_region(function(r) fread(file.path(qc_dir(r), "qc", "genomic_context.tsv")))
gdist    <- read_by_region(function(r) fread(file.path(qc_dir(r), "qc", "distance_to_nearest_gene.tsv")))

for (d in list(vmr, cutoffs, tqc, coverage, turnover, gctx, gdist)) d[, region := as_region(region)]

SCRIPT <- "11_integrated_manuscript_outputs/_h/01_figure1_catalog.R"
runs_used <- c(vapply(regions, CATALOG_RUN, ""), vapply(regions, QC_RUN, ""))

## --------------------------------------------- A. donors and CpGs per region
##
## Autosomes only: sex chromosomes are prepared but held out of the primary
## catalog (V4), so including them would inflate the assayed-CpG count.
design <- tqc[is_primary_chrom == TRUE,
              .(n_donors = max(n_donors), n_cpgs = sum(n_cpgs),
                n_vmrs = NA_integer_), by = region]
design[vmr[, .N, by = region], n_vmrs := i.N, on = "region"]

design_long <- melt(design, id.vars = "region",
                    measure.vars = c("n_donors", "n_cpgs", "n_vmrs"),
                    variable.name = "metric", value.name = "value")
design_long[, metric := factor(metric, levels = c("n_donors", "n_cpgs", "n_vmrs"),
                               labels = c("Donors", "CpGs assayed", "VMRs called"))]
design_long[, lab := fifelse(
    metric == "Donors", label_comma(accuracy = 1)(value),
    label_number(scale_cut = cut_short_scale(), accuracy = 0.1)(value))]

pA <- ggplot(design_long, aes(region, value, fill = region)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = lab), vjust = -0.35, size = 2.5, colour = "black") +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = REGION_COLORS, guide = "none") +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = NULL) +
    BASE_THEME + NO_TITLES +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

## ------------------------------------------- B. VMRs per chromosome, overlaid
chr_levels <- paste0("chr", chrom_order(include_sex = FALSE))
cutoffs[, chr := factor(chr, levels = chr_levels)]

pB <- ggplot(cutoffs[!is.na(chr)], aes(chr, n_vmrs, colour = region, group = region)) +
    geom_line(linewidth = 0.45, alpha = 0.85) +
    geom_point(size = 1.15) +
    scale_colour_manual(values = REGION_COLORS, name = NULL) +
    scale_x_discrete(labels = function(x) sub("^chr", "", x)) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
    labs(x = "Chromosome", y = "VMRs called") +
    BASE_THEME + NO_TITLES +
    ## ggplot2 3.5 split the inside-legend position out of `legend.position`;
    ## the bare c(x, y) form is soft-deprecated in 4.0 (this env runs 4.0.1).
    theme(legend.position = "inside",
          legend.position.inside = c(0.86, 0.82),
          legend.background = element_blank(),
          legend.key.height = grid::unit(0.34, "cm"))

## ------------------------------------------ C. VMR width and CpGs per VMR
vmr[, width := end - start + 1L]
shape <- rbind(
    vmr[, .(region, metric = "VMR width (bp)", value = as.numeric(width))],
    vmr[, .(region, metric = "CpGs per VMR", value = as.numeric(n))])
shape[, metric := factor(metric, levels = c("VMR width (bp)", "CpGs per VMR"))]

pC <- ggplot(shape, aes(value, colour = region)) +
    geom_density(linewidth = 0.5) +
    facet_wrap(~ metric, scales = "free", nrow = 1) +
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    labs(x = NULL, y = "Density") +
    BASE_THEME + NO_TITLES

## ------------------------------------- D. coverage outside array-accessible CpGs
##
## The comparator platform is named on the axis, not only in the caption: the
## fraction is meaningless without knowing which array it is relative to.
cov_p <- coverage[array_platform == platform]
if (nrow(cov_p) == 0) {
    stop("No array coverage rows for platform ", platform,
         ". Run 01_vmr_catalog/_h/04b_rerun_array_coverage.R first.")
}

cov_long <- melt(cov_p, id.vars = "region",
                 measure.vars = c("frac_vmr_cpgs_off_array",
                                  "frac_vmrs_invisible_to_array"),
                 variable.name = "metric", value.name = "frac")
cov_long[, metric := factor(
    metric,
    levels = c("frac_vmr_cpgs_off_array", "frac_vmrs_invisible_to_array"),
    labels = c("CpGs", "VMRs"))]

pD <- ggplot(cov_long, aes(region, frac, fill = region)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = percent(frac, accuracy = 0.1)),
              vjust = -0.35, size = 2.5, colour = "black") +
    facet_wrap(~ metric, nrow = 1) +
    scale_fill_manual(values = REGION_COLORS, guide = "none") +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1),
                       expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL,
         y = paste0("Fraction outside\nIllumina ", platform, " coverage")) +
    BASE_THEME + NO_TITLES +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

## -------------------------------- supplement: catalog turnover (audit only)
##
## The legacy catalog is a comparison baseline only (AGENTS.md 8); it is
## invalid for scientific use because of the V1 donor-row misalignment. This
## shows how much the repair moved, not that either set is correct -- which is
## why it is a supplementary figure and not a main panel.
turn <- turnover[, .(region,
                     Retained = n_v2_overlapping_legacy,
                     `Novel in v2` = n_v2_novel,
                     `Lost from legacy` = n_legacy_lost)]
turn_long <- melt(turn, id.vars = "region", variable.name = "class",
                  value.name = "n")
turn_long[, class := factor(class, levels = c("Retained", "Novel in v2",
                                              "Lost from legacy"))]

pTurn <- ggplot(turn_long, aes(region, n, fill = class)) +
    geom_col(width = 0.8, position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = c(Retained = PAL_CHARCOAL,
                                 `Novel in v2` = PAL_RUST,
                                 `Lost from legacy` = PAL_NULL), name = NULL) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       expand = expansion(mult = c(0, 0.10))) +
    labs(x = NULL, y = "VMRs") +
    BASE_THEME + NO_TITLES +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0))

## ------------------------------------------- F. genomic context of the catalog
##
## Descriptive only: compartments are assigned by priority so they partition
## the catalog. No enrichment is claimed here (that is Module 04's remit, and
## the enrichment model is a locked PI decision -- AGENTS.md 12).
CTX_LEVELS <- c("Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Intergenic")
gctx[, genomic_context := factor(genomic_context, levels = rev(CTX_LEVELS))]

pF <- ggplot(gctx, aes(frac_vmrs, genomic_context, fill = region)) +
    geom_col(width = 0.7, position = position_dodge(width = 0.75)) +
    scale_fill_manual(values = REGION_COLORS, name = NULL) +
    scale_x_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(x = "VMRs", y = NULL) +
    BASE_THEME + NO_TITLES +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0))

## ------------------------------------------- G. distance to the nearest gene
##
## Zero-distance VMRs (those inside a gene) cannot be shown on a log axis, so
## they are reported as an explicit annotation rather than dropped silently.
frac_in_gene <- gdist[, .(f = mean(distance_to_nearest_gene == 0)), by = region]
pG <- ggplot(gdist[distance_to_nearest_gene > 0],
             aes(distance_to_nearest_gene, colour = region)) +
    geom_density(linewidth = 0.5) +
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    annotate("text", x = 1, y = Inf, hjust = 0, vjust = 1.6, size = 2.4,
             colour = "grey35",
             label = sprintf("%s%% of VMRs lie within a gene body",
                             paste(sprintf("%.0f", 100 * frac_in_gene$f),
                                   collapse = "/"))) +
    labs(x = "Distance to nearest gene (bp)", y = "Density") +
    BASE_THEME + NO_TITLES

## ------------------------------------------------------------------ assemble
##
## Four rows, <= FIG_HEIGHT_MAX. save_figure() enforces the ceiling; it is not
## left to a visual check.
row3 <- pC + pD + plot_layout(widths = c(1, 0.85))
row4 <- pF + pG + plot_layout(widths = c(1, 1))
figure <- (pA / pB / row3 / row4) +
    plot_layout(heights = c(1.0, 1.05, 1.15, 1.05)) +
    fig_tags(theme = theme(plot.margin = margin(2, 2, 2, 2))) &
    TAG_THEME

suffix <- if (platform == "450K") "" else paste0("_", tolower(platform))
arm <- if (cohort == "AA") "" else paste0("_", cohort)
STEM <- paste0("figure1_vmr_catalog", arm, suffix)
save_figure(figure, STEM, width = FIG_WIDTH_FULL, height = 9.2, fig_dir = fig_dir)

## The turnover audit, as its own supplementary figure.
##
## Written on the 450K pass only. Turnover compares v2 VMR calls against the
## legacy ones and has nothing to do with which array is the comparator, so the
## EPIC pass was writing a byte-identical second copy under a different name --
## four turnover figures where there are two results.
TURN_STEM <- paste0("figureS_catalog_turnover", arm)
if (platform == "450K") {
    save_figure(pTurn + fig_tags() & TAG_THEME, TURN_STEM,
                width = FIG_WIDTH_THREEQ, height = 3.2, fig_dir = fig_dir)
}

## ---------------------------------------------------------- source data
##
## Panel names follow the RENDERED tag, read off the assembly order above, not
## the R variable holding the plot. The previous version named gctx "panelF"
## and gdist "panelG" while patchwork tagged them e and f, so every panel's
## source data documented the wrong panel.
PANELS <- c(a = "pA", b = "pB", c = "pC", d = "pD", e = "pF", f = "pG")
tag_for <- function(var) {
    hit <- names(PANELS)[PANELS == var]
    if (length(hit) != 1) stop("No rendered tag for plot object ", var)
    hit
}

sd <- function(dt, var, tbl, filt) {
    write_source_data(dt, paste0(STEM, "_panel_", tag_for(var)),
                      runs_used, tbl, SCRIPT, filt, data_dir)
}
sd(design_long, "pA", "qc/technical_qc.tsv + vmr/vmr_catalog.tsv",
   "is_primary_chrom == TRUE (autosomes; sex chromosomes excluded per V4)")
sd(cutoffs[!is.na(chr), .(region, chr, n_vmrs, sd_cutoff, n_cpgs_tested)],
   "pB", "vmr/sd_cutoffs.tsv", "autosomes only")
sd(shape[, .(n = .N, median = median(value), q25 = quantile(value, .25),
             q75 = quantile(value, .75)), by = .(region, metric)],
   "pC", "vmr/vmr_catalog.tsv", "all called VMRs; summary of plotted density")
sd(cov_long, "pD", "qc/array_coverage.tsv",
   paste0("array_platform == '", platform, "'"))
sd(gctx, "pF", "qc/genomic_context.tsv",
   "all called VMRs; compartments assigned by priority, so they partition the catalog")
sd(gdist[, .(n = .N, frac_within_gene = mean(distance_to_nearest_gene == 0),
             median_nonzero = as.numeric(median(distance_to_nearest_gene[distance_to_nearest_gene > 0]))),
         by = region],
   "pG", "qc/distance_to_nearest_gene.tsv",
   "all called VMRs; density plotted for distance > 0 only (log axis)")

## The turnover audit ships as its own figure, so its source data is named for
## that figure rather than for Figure 1.
if (platform == "450K") {
    write_source_data(turn_long, paste0(TURN_STEM, "_panel_a"), runs_used,
                      "qc/vmr_turnover.tsv", SCRIPT,
                      "autosomes; legacy is a comparison baseline only (AGENTS.md 8)",
                      data_dir)
}

message("[done] Figure 1 written to ", fig_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
