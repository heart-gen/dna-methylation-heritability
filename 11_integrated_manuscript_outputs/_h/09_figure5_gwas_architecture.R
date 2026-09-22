#### 11 / Figure 5: GWAS loci and the local genetic-control axis ####
##
## The axis depletion is TRAIT-GENERAL. Module 09 stage 17 ran the identical
## locus -> VMR -> axis contrast over the harmonized GWAS collection under one
## lead-SNP rule, and schizophrenia sits at the 16th-33rd percentile of the 63
## distribution traits depending on region, with psychiatric traits
## indistinguishable as a category (Wilcoxon p 0.11-0.96).
##
## So the main figure is the DISTRIBUTION, not schizophrenia. AGENTS.md 7.8
## rule 1 and 11 both require the depletion be written as a property of
## trait-associated loci with schizophrenia as a typical example, and a figure
## that gave SCZ its own panels would contradict the text it illustrates. The
## schizophrenia locus evidence -- prioritized loci, meQTL support, coupling,
## colocalization -- is figureS_schizophrenia_application.
##
## Module 09's decision 2 is scz_application_retention = RETAIN_MAIN_TEXT, and
## that is honoured: schizophrenia appears in panels a, b and c as a marked
## example among the other traits, never as a privileged encoding.
##
## Panels, in rendered order
##   a  every trait's axis estimate, grouped by GWAS category, three regions,
##      with schizophrenia marked in place
##   b  schizophrenia against the null distribution it belongs to
##   c  psychiatric vs other traits -- indistinguishable as a category
##   d  the mechanism: what a trait's depletion actually tracks
##
## WHAT THIS FIGURE MAY NOT SAY
##  - never "schizophrenia-specific" (AGENTS.md 7.8 rule 1);
##  - never a repeat result. NO trait reaches q < 0.05 for LINE/L1 enrichment
##    in any region or arm; the LINE/L1 axis-link correlation is inconsistent
##    in sign across regions and collapses under the high-mappability
##    restriction, so panel d shows BOTH arms rather than the headline one;
##  - the GWAS collection is European or European-dominated while the cohort is
##    admixed African American. The limitation attaches to the LOCUS
##    DEFINITION, not to the axis, which is a within-cohort rank;
##  - stages 17/18 QUALIFY Module 09; they never validate it;
##  - the extended MHC is excluded for every trait, which removes the immune
##    category's strongest locus;
##  - caudate is batch-confounded (AGENTS.md 8.1).
##
## Usage:
##   Rscript 09_figure5_gwas_architecture.R --cohort AA --run-id fig-all-YYYYMMDD
##   Rscript 09_figure5_gwas_architecture.R --cohort AA --run-id smoke \
##       --out-dir /path/to/draft

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(V2_ROOT, "11_integrated_manuscript_outputs", "_h",
                 "00_figure_theme.R"))

suppressPackageStartupMessages({
    library(patchwork)
    library(scales)
})

opts <- parse_v2_args(require = c("cohort", "run_id"))
cohort  <- opts$cohort
regions <- load_config("cohorts")$regions

module_root <- file.path(V2_ROOT, "11_integrated_manuscript_outputs")
run_dir  <- if (!is.null(opts$out_dir)) opts$out_dir else
    file.path(module_root, "_m", "runs", opts$run_id)
fig_dir  <- file.path(run_dir, "figures")
data_dir <- file.path(run_dir, "source_data")

SCRIPT <- "11_integrated_manuscript_outputs/_h/09_figure5_gwas_architecture.R"

## ------------------------------------------------------------- upstream
##
## Stages 17 and 18 are module-level: they read the accepted per-region runs
## and write _m/combined/. Resolving the per-region runs through the gate is
## still what certifies the tables were built on accepted input.
SCZ <- vapply(regions, function(r)
    require_accepted_upstream("09_schizophrenia_risk_application", cohort, r)$run_id,
    character(1))
comb <- file.path(V2_ROOT, "09_schizophrenia_risk_application", "_m", "combined")

need <- function(f) {
    p <- file.path(comb, f)
    if (!file.exists(p)) {
        stop("Missing ", f, ". Stages 17/18 must be re-run WITHOUT ",
             "--allow-unlocked so their output is citable: ",
             "sbatch 09_schizophrenia_risk_application/_h/step_9_negative_controls.sh")
    }
    fread(p)
}
traits  <- need(sprintf("scz-negative-control-traits-%s.tsv", cohort))
summ    <- need(sprintf("scz-negative-control-summary-%s.tsv", cohort))
bycat   <- need(sprintf("scz-negative-control-by-category-%s.tsv", cohort))
axlink  <- need(sprintf("scz-locus-architecture-axis-link-%s.tsv", cohort))
arch    <- need(sprintf("scz-locus-architecture-%s.tsv", cohort))

## -------------------------------------------------------------- guards
##
## An -UNACCEPTED file may never reach a figure (AGENTS.md 6), and the files
## above are the unsuffixed ones by construction -- but the flag is checked
## rather than assumed.
for (nm in c("traits", "summ", "axlink", "arch")) {
    d <- get(nm)
    if ("built_with_unaccepted_runs" %in% names(d) &&
        any(d$built_with_unaccepted_runs %in% c(TRUE, "TRUE"))) {
        stop(nm, " carries built_with_unaccepted_runs = TRUE and may not be ",
             "cited (AGENTS.md 6).")
    }
}

## The headline is trait-GENERAL. If schizophrenia ever stopped being typical,
## this figure's whole framing would be wrong, so the claim is checked.
PCT <- summ[, range(scz_percentile_in_distribution, na.rm = TRUE)]
WIL <- summ[, min(psychiatric_vs_other_wilcoxon_p, na.rm = TRUE)]
if (WIL < 0.05) {
    stop("Psychiatric traits now separate from the rest (min Wilcoxon p = ",
         signif(WIL, 3), "). Figure 5 is built on their being ",
         "indistinguishable; re-read stage 17 before rebuilding.")
}
if (PCT[1] < 0.02) {
    stop("Schizophrenia is now extreme in the trait distribution (percentile ",
         signif(PCT[1], 3), "), so presenting it as a typical example would ",
         "misstate the result.")
}

## NO trait may reach q < 0.05 for LINE/L1 ENRICHMENT, in any region or arm.
## This is the per-trait annotation test, NOT the axis-link Spearman -- the
## axis-link correlation does have nominally significant LINE/L1 rows in the
## all_vmrs arm, in inconsistent directions, which is exactly why panel d shows
## the high-mappability arm beside it.
l1 <- arch[annotation %like% "line_l1" & !(skipped %in% c(TRUE, "TRUE"))]
if (any(l1$q < 0.05, na.rm = TRUE)) {
    stop("A trait now reaches q < 0.05 for LINE/L1 enrichment (", 
         sum(l1$q < 0.05, na.rm = TRUE), " of ", nrow(l1), " tests). ",
         "Panel d is captioned on there being none; re-read stage 18.")
}
message("[guard] SCZ percentile ", paste(round(PCT, 3), collapse = "-"),
        "; psychiatric-vs-other Wilcoxon p >= ", signif(WIL, 2),
        "; LINE/L1 enrichment 0/", nrow(l1), " tests at q < 0.05")

PRIMARY_MODEL <- "A_technical"
MODEL_LABELS <- c(A_technical = "Technical covariates",
                  B_technical_context = "+ genomic context")

## ------------------------------- a. every trait, grouped by GWAS category
dist <- traits[in_distribution == TRUE & model == PRIMARY_MODEL &
               is.finite(estimate)]
dist[, region := as_region(region)]
## Categories ordered by how depleted they are, pooling regions, so the reader
## reads the panel top-to-bottom as a gradient.
ord <- dist[, .(m = median(estimate)), by = category][order(m)]
dist[, category := factor(category, levels = rev(ord$category))]

scz <- summ[model == PRIMARY_MODEL, .(region, scz_estimate, distribution_median,
                                      distribution_q25, distribution_q75,
                                      scz_percentile_in_distribution,
                                      n_traits_in_distribution, scz_p)]
scz[, region := as_region(region)]

pA <- ggplot(dist, aes(estimate, category)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    geom_vline(data = scz, aes(xintercept = scz_estimate),
               colour = PAL_CHARCOAL, linetype = "dashed", linewidth = 0.4) +
    geom_point(aes(colour = region), size = 0.85, alpha = 0.75) +
    facet_wrap(~ region, nrow = 1) +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    labs(x = "Axis estimate: trait-linked VMRs vs background (score SD)",
         y = NULL,
         caption = paste0(
             "One point per trait (n = ", scz$n_traits_in_distribution[1],
             " traits with 10 or more leads). Dashed line, schizophrenia.\n",
             "GWAS loci defined from European or European-dominated summary",
             " statistics.\nExtended MHC excluded for every trait."))  +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

## ----------------------------- b. schizophrenia against its own distribution
pB <- ggplot(scz, aes(y = region)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    errorbar_h(aes(xmin = distribution_q25, xmax = distribution_q75),
               colour = PAL_NULL, linewidth = 1.8, alpha = 0.8) +
    geom_point(aes(x = distribution_median), shape = 124, size = 2.6,
               colour = PAL_CHARCOAL) +
    geom_point(aes(x = scz_estimate, colour = region), size = 2.1) +
    ## Anchored to the panel's right edge, not to the point: at the left end of
    ## the scale the label ran off the panel entirely.
    geom_text(aes(label = paste0(percent(scz_percentile_in_distribution,
                                         accuracy = 1), " of traits lower")),
              x = Inf, hjust = 1.05, vjust = -1.1, size = 2.2,
              colour = "grey35") +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_y_discrete(limits = rev, expand = expansion(add = c(0.5, 0.8))) +
    labs(x = "Axis estimate (score SD)", y = NULL,
         caption = paste("Grey bar, trait-distribution IQR;\ntick, its",
                         "median; point, schizophrenia.")) +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

## --------------------------- c. psychiatric traits are not a separate class
psy <- traits[in_distribution == TRUE & model == PRIMARY_MODEL &
              is.finite(estimate)]
## Short labels: "Psychiatric-neurologic" and "All other traits" overlapped at
## this panel's width and squeezed the facet strips into "ppocamp". The full
## grouping is in the caption and the source data.
psy[, grp := fifelse(category == "Psychiatric-neurologic",
                     "Psychiatric", "Other")]
psy[, grp := factor(grp, levels = c("Psychiatric", "Other"))]
psy[, region := as_region(region)]
wil <- summ[model == PRIMARY_MODEL, .(region, psychiatric_vs_other_wilcoxon_p)]
wil[, region := as_region(region)]
wil[, lab := paste0("P = ", signif(psychiatric_vs_other_wilcoxon_p, 2))]

## Region on x, group as fill, and NO facets. Faceting by region gave each
## panel about 1.2 in for two category labels and a strip, so the labels
## collided and "Hippocampus" was clipped to "ppocampu". One panel also puts
## the two groups side by side within a region, which is the comparison.
pC <- ggplot(psy, aes(region, estimate, fill = grp)) +
    geom_hline(yintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    geom_boxplot(width = 0.62, outlier.size = 0.4, linewidth = 0.3,
                 position = position_dodge(width = 0.72)) +
    geom_text(data = wil, aes(x = region, y = Inf, label = lab),
              inherit.aes = FALSE, vjust = 1.5, size = 2.3, colour = "grey35") +
    scale_fill_manual(values = c(Psychiatric = PAL_TAN, Other = PAL_NULL),
                      name = NULL) +
    ## Headroom so the P label clears the topmost outlier.
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
    labs(x = NULL, y = "Axis estimate (score SD)",
         caption = paste("Psychiatric = psychiatric-neurologic;",
                         "Other = all remaining.\nTwo-sided Wilcoxon.")) +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

## ------------------------------------- d. what a trait's depletion tracks
ANNOT_LABELS <- c(accessible_any = "Accessible", h3k27ac_any = "H3K27ac",
                  bivalent_any = "Bivalent", h3k27me3_any = "H3K27me3",
                  h3k9me3_any = "H3K9me3", quiescent_any = "Quiescent",
                  line_l1_any = "LINE/L1")
ARM_LABELS <- c(all_vmrs = "All VMRs", high_mappability = "High mappability")

ax <- axlink[adjustment == "technical" & annotation %in% names(ANNOT_LABELS)]
ax[, region := as_region(region)]
ax[, alab := factor(ANNOT_LABELS[annotation], levels = rev(unname(ANNOT_LABELS)))]
ax[, armf := factor(ARM_LABELS[arm], levels = unname(ARM_LABELS))]
ax[, stars := sig_stars(q)]

pD <- ggplot(ax, aes(spearman_rho, alab, colour = region)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    geom_point(size = 1.5, position = position_dodge(width = 0.66)) +
    geom_text(aes(label = stars), position = position_dodge(width = 0.66),
              hjust = -0.4, vjust = 0.75, size = 2.5, show.legend = FALSE) +
    facet_wrap(~ armf, nrow = 1) +
    scale_colour_manual(values = REGION_COLORS, name = NULL) +
    scale_x_continuous(limits = c(-0.75, 0.75), breaks = seq(-0.5, 0.5, 0.5)) +
    labs(x = "Spearman: trait's axis depletion vs its loci's annotation enrichment",
         y = NULL,
         caption = paste0(
             "No trait reaches q < 0.05 for LINE/L1 ENRICHMENT in any region",
             " or arm (0 of ", nrow(l1), " tests).\nThe correlation shown is",
             " inconsistent in sign and does not survive the high-mappability",
             " arm.\nOverlap is overlap: no activity or retrotransposition",
             " claim follows."))  +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

## ------------------------------------------------------------------ assemble
arm_sfx <- if (cohort == "AA") "" else paste0("_", cohort)
STEM <- paste0("figure5_gwas_architecture_axis", arm_sfx)

figure <- (pA / ((pB | pC) + plot_layout(widths = c(0.8, 1))) / pD) +
    plot_layout(heights = c(1.25, 0.85, 1.05)) +
    fig_tags() & TAG_THEME

save_figure(figure, STEM, width = FIG_WIDTH_FULL, height = 9.0,
            fig_dir = fig_dir)

## ---------------------------------------------------------- source data
runs_used <- unname(SCZ)
sd <- function(dt, nm, tbl, filt) {
    write_source_data(dt, paste0(STEM, "_", nm), runs_used, tbl, SCRIPT,
                      filt, data_dir)
}
ANC <- paste("GWAS loci defined from European or European-dominated summary",
             "statistics; the methylation cohort is admixed African American.",
             "The limitation attaches to the locus definition, not to the axis,",
             "which is a within-cohort rank. Extended MHC excluded for every trait.")

sd(dist[, .(region, model, tag, category, role, n_leads, n_linked, n_background,
            estimate, se, p, q)],
   "panel_a", "_m/combined/scz-negative-control-traits-{cohort}.tsv",
   paste("in_distribution == TRUE; model ==", PRIMARY_MODEL, ".", ANC))
sd(scz, "panel_b", "_m/combined/scz-negative-control-summary-{cohort}.tsv",
   paste("model ==", PRIMARY_MODEL,
         "; schizophrenia read against the trait distribution it belongs to.",
         "Trait-general: SCZ is a typical example, never schizophrenia-specific."))
sd(wil, "panel_c", "_m/combined/scz-negative-control-summary-{cohort}.tsv",
   paste("model ==", PRIMARY_MODEL,
         "; psychiatric vs all other traits, two-sided Wilcoxon rank-sum"))
sd(ax[, .(region, annotation, adjustment, arm, n_traits, spearman_rho,
          spearman_p, q, median_annotation_estimate, interpretation)],
   "panel_d", "_m/combined/scz-locus-architecture-axis-link-{cohort}.tsv",
   paste("adjustment == 'technical'; both arms shown.",
         "NO trait reaches q < 0.05 for LINE/L1 enrichment in any region or arm;",
         "the axis-link LINE/L1 correlation is inconsistent in sign and collapses",
         "under the high-mappability restriction (AGENTS.md 7.8)."))
sd(bycat, "category_medians",
   "_m/combined/scz-negative-control-by-category-{cohort}.tsv",
   "category medians for both models; broad genomic context moves the estimate 1-3%")
sd(summ, "trait_distribution_summary",
   "_m/combined/scz-negative-control-summary-{cohort}.tsv",
   "both models; stages 17/18 QUALIFY Module 09 and change neither of its decisions")

message("[done] Figure 5 written to ", fig_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
