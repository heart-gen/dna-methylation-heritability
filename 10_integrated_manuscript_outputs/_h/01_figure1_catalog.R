#### 10 / Figure 1: cohort, corrected VMR catalog, and off-array coverage ####
##
## AGENTS.md 11 assigns Figure 1 the cohort, the corrected VMR catalog, and WGBS
## coverage outside array-accessible CpGs. AGENTS.md 7.9 makes this module the
## only place manuscript figures may be assembled.
##
## Panels
##   A  donors and CpGs assayed per region
##   B  VMRs per chromosome, three regions
##   C  VMR width and CpGs per VMR
##   D  WGBS coverage outside array-accessible CpGs   <- the 2.2 contribution
##   E  genomic compartment the VMRs occupy (descriptive, not enrichment)
##   F  distance to the nearest gene
##   G  catalog turnover against the invalid legacy calls (the V1 repair)
##
## Usage:
##   Rscript 01_figure1_catalog.R --cohort AA --run-id fig-all-20260826
##   Rscript 01_figure1_catalog.R --cohort AA --run-id ID --platform EPIC

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(V2_ROOT, "10_integrated_manuscript_outputs", "_h",
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

module_root <- file.path(V2_ROOT, "10_integrated_manuscript_outputs")
run_dir  <- file.path(module_root, "_m", "runs", opts$run_id)
fig_dir  <- file.path(run_dir, "figures")
data_dir <- file.path(run_dir, "source_data")

## The catalog run holds the VMR calls; the qc rerun holds array coverage,
## which the sealed catalog runs predate.
CATALOG_RUN <- function(r) paste0("vmrcat-", cohort, "-", r, "-20260816")
QC_RUN      <- function(r) paste0("vmrcatqc-", cohort, "-", r, "-20260826-a")

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

SCRIPT <- "10_integrated_manuscript_outputs/_h/01_figure1_catalog.R"
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
    theme(legend.position = c(0.86, 0.82),
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

## ---------------------------------------------------- E. catalog turnover
##
## The legacy catalog is a comparison baseline only (AGENTS.md 8); it is
## invalid for scientific use because of the V1 donor-row misalignment. This
## panel shows how much the repair moved, not that either set is correct.
turn <- turnover[, .(region,
                     Retained = n_v2_overlapping_legacy,
                     `Novel in v2` = n_v2_novel,
                     `Lost from legacy` = n_legacy_lost)]
turn_long <- melt(turn, id.vars = "region", variable.name = "class",
                  value.name = "n")
turn_long[, class := factor(class, levels = c("Retained", "Novel in v2",
                                              "Lost from legacy"))]

pE <- ggplot(turn_long, aes(region, n, fill = class)) +
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
row3 <- pC + pD + plot_layout(widths = c(1, 0.85))
row4 <- pF + pG + plot_layout(widths = c(1, 1))
figure <- (pA / pB / row3 / row4 / pE) +
    plot_layout(heights = c(1.0, 1.05, 1.15, 1.05, 0.9)) +
    plot_annotation(tag_levels = "A",
                    theme = theme(plot.margin = margin(2, 2, 2, 2))) &
    theme(plot.tag = element_text(face = "bold", size = 11))

suffix <- if (platform == "450K") "" else paste0("_", tolower(platform))
arm <- if (cohort == "AA") "" else paste0("_", cohort)
save_figure(figure, paste0("figure1_vmr_catalog", arm, suffix),
            width = FIG_WIDTH_FULL, height = 11.4, fig_dir = fig_dir)

## ---------------------------------------------------------- source data
sd <- function(dt, nm, tbl, filt) {
    write_source_data(dt, paste0("figure1_vmr_catalog", arm, suffix, "_", nm),
                      runs_used, tbl, SCRIPT, filt, data_dir)
}
sd(design_long, "panelA", "qc/technical_qc.tsv + vmr/vmr_catalog.tsv",
   "is_primary_chrom == TRUE (autosomes; sex chromosomes excluded per V4)")
sd(cutoffs[!is.na(chr), .(region, chr, n_vmrs, sd_cutoff, n_cpgs_tested)],
   "panelB", "vmr/sd_cutoffs.tsv", "autosomes only")
sd(shape[, .(n = .N, median = median(value), q25 = quantile(value, .25),
             q75 = quantile(value, .75)), by = .(region, metric)],
   "panelC", "vmr/vmr_catalog.tsv", "all called VMRs; summary of plotted density")
sd(cov_long, "panelD", "qc/array_coverage.tsv",
   paste0("array_platform == '", platform, "'"))
sd(gctx, "panelF", "qc/genomic_context.tsv",
   "all called VMRs; compartments assigned by priority, so they partition the catalog")
sd(gdist[, .(n = .N, frac_within_gene = mean(distance_to_nearest_gene == 0),
             median_nonzero = as.numeric(median(distance_to_nearest_gene[distance_to_nearest_gene > 0]))),
         by = region],
   "panelG", "qc/distance_to_nearest_gene.tsv",
   "all called VMRs; density plotted for distance > 0 only (log axis)")
sd(turn_long, "panelE", "qc/vmr_turnover.tsv",
   "autosomes; legacy is a comparison baseline only (AGENTS.md 8)")

message("[done] Figure 1 written to ", fig_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
