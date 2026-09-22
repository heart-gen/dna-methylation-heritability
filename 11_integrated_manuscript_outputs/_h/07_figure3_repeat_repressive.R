#### 11 / Figure 3: repeat-rich and repressive chromatin architecture ####
##
## AGENTS.md 11 assigns Figure 3 "LINE/L1, H3K9me3, quiescent chromatin,
## mappability, and cell sensitivities". Module 04 is the owning analysis
## (AGENTS.md 7.4) and its accepted runs are rra-AA-{region}-20260906.
##
## Panels, in rendered order
##   a  the BH family: H3K9me3, LINE/L1, quiescent, three regions
##   b  complementary contrasts (accessible, H3K27ac) and the H3K27me3
##      specificity control -- the compartments that run the OTHER way
##   c  the five locked analysis sets, so a reader sees the sensitivities
##      rather than being told they passed
##   d  the continuous gradient behind panel a, across score deciles
##
## WHAT THIS FIGURE MAY AND MAY NOT SAY
## The permitted claim per outcome is read from the run at build time and
## asserted, not restated from memory. As accepted on 2026-09-08:
##   quiescent_frac  shared across all three regions (3/3)
##   h3k9me3_frac    BELOW the gate (2/3); suggestive in dlpfc + hippocampus
##                   only, and never described as shared
##   line_l1_frac    supported in both ELIGIBLE regions (dlpfc, hippocampus).
##                   Caudate is set aside as technically confounded -- it is
##                   NOT a null result and NOT a failure to replicate, and it
##                   is excluded from the 2/2 denominator. It renders in
##                   PAL_NULL with an explicit "set aside" mark.
## Caudate is batch-confounded throughout (AGENTS.md 8.1), so no caudate-vs-
## other-region difference anywhere in this figure is attributable to region.
## Overlap is overlap: no activity, expression or retrotransposition claim
## follows from any panel here (AGENTS.md 2.3, 7.4).
##
## `outcome_role` is carried, never flattened. The BH family is three outcomes;
## accessible/H3K27ac are complementary contrasts, H3K27me3 is a specificity
## control and bivalent is descriptive. They do not share an error rate and
## must not share a panel with the BH family as if they did.
##
## Usage:
##   Rscript 07_figure3_repeat_repressive.R --cohort AA --run-id fig-all-YYYYMMDD
##   Rscript 07_figure3_repeat_repressive.R --cohort AA --run-id smoke \
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

SCRIPT <- "11_integrated_manuscript_outputs/_h/07_figure3_repeat_repressive.R"
PREDICTOR <- "local_snp_contribution_score_z"

## ------------------------------------------------------------- upstream
RRA <- vapply(regions, function(r)
    require_accepted_upstream("04_repeat_repressive_architecture", cohort, r)$run_id,
    character(1))
rra_dir <- function(r) file.path(V2_ROOT, "04_repeat_repressive_architecture",
                                 "_m", "runs", RRA[[r]], "results")

## The claim table and the cross-region association table live on the GATE HOST
## (the caudate run), not on each region's own run. Reading them per region
## would silently give three copies of the caudate view.
GATE_HOST <- "caudate"
claims <- fread(file.path(rra_dir(GATE_HOST), "interpretation-claims.tsv"))
assoc  <- fread(file.path(rra_dir(GATE_HOST), "association-results-all-regions.tsv"))

## -------------------------------------------------------------- guards
##
## A figure must not be able to assert more than the accepted decision allows,
## so the contract is checked here rather than trusted to the panel code.
banned <- intersect(c("h2_en_calibrated", "r_squared_cv", "h2_unscaled"),
                    names(assoc))
if (length(banned) > 0) {
    stop("Retired quantity in the Module 04 contract: ",
         paste(banned, collapse = ", "), " (AGENTS.md 3).")
}
if (!setequal(unique(assoc$region), regions)) {
    stop("association-results-all-regions.tsv covers ",
         paste(unique(assoc$region), collapse = "/"), ", not all three regions.")
}

## Ordered strongest-claim-first, and this order is used in EVERY panel: a's
## y-axis, c's facets and d's facets. Panels that rank the same three outcomes
## differently make a reader re-learn the figure at each row.
BH_FAMILY <- c("quiescent_frac", "h3k9me3_frac", "line_l1_frac")
stopifnot(setequal(claims$outcome, BH_FAMILY))
stopifnot(setequal(
    unique(assoc[outcome %in% BH_FAMILY, outcome_role]), "bh_family"))

## LINE/L1 in caudate is SET ASIDE, not null. If that ever stops being true the
## panel's encoding is wrong and the figure must stop rather than quietly
## promote caudate into the claim.
l1 <- claims[outcome == "line_l1_frac"]
SET_ASIDE <- trimws(strsplit(l1$regions_excluded, ",")[[1]])
SET_ASIDE <- SET_ASIDE[nzchar(SET_ASIDE)]
if (!identical(SET_ASIDE, "caudate")) {
    stop("Module 04 no longer sets caudate aside for LINE/L1 (regions_excluded = '",
         l1$regions_excluded, "'). Figure 3 encodes that exclusion explicitly; ",
         "re-read the accepted claim before rebuilding.")
}
stopifnot(l1$regions_eligible == 2L, l1$regions_surviving == 2L)

## H3K9me3 is below the shared gate. Never rendered as a 3/3 result.
h9 <- claims[outcome == "h3k9me3_frac"]
stopifnot(h9$regions_surviving < h9$regions_required)
qu <- claims[outcome == "quiescent_frac"]
stopifnot(qu$regions_surviving == qu$regions_required)

message("[guard] claims as accepted: quiescent ", qu$regions_surviving, "/",
        qu$regions_required, "; h3k9me3 ", h9$regions_surviving, "/",
        h9$regions_required, " (below gate); line_l1 ", l1$regions_surviving,
        "/", l1$regions_eligible, ", caudate set aside")

## ------------------------------------------------------------- labelling
OUTCOME_LABELS <- c(
    quiescent_frac  = "Quiescent",
    h3k9me3_frac    = "H3K9me3",
    line_l1_frac    = "LINE/L1",
    accessible_frac = "Accessible",
    h3k27ac_frac    = "H3K27ac",
    h3k27me3_frac   = "H3K27me3",
    bivalent_frac   = "Bivalent")

## Kept short for the same reason as OUTCOME_LABELS: patchwork sizes the whole
## column's left margin from the longest axis label in it. The full analysis-set
## names are in the source data.
SET_LABELS <- c(
    primary                 = "Primary",
    high_mappability        = "High mappability",
    exclude_segdups         = "No segdups",
    adjust_cell_composition = "Cell adjusted",
    low_cell_composition    = "Low cell content")

prep <- function(dt) {
    d <- copy(dt)
    d[, region := as_region(region)]
    ## `label` orders a DISCRETE Y AXIS, where the first level sits at the
    ## bottom, so it is reversed. `flab` orders FACETS left to right, so it is
    ## not. Both come from the same OUTCOME_LABELS order.
    d[, label := factor(OUTCOME_LABELS[outcome],
                        levels = rev(unname(OUTCOME_LABELS)))]
    d[, flab := factor(OUTCOME_LABELS[outcome], levels = unname(OUTCOME_LABELS))]
    d[, `:=`(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)]
    d[]
}

## -------------------------------------------- a. the BH family, three regions
prim <- prep(assoc[analysis_set == "primary" & predictor == PREDICTOR &
                   outcome %in% BH_FAMILY])
## Grey is reserved for "set aside from the claim", and the mark says so. A
## reader must not be able to mistake it for a null.
prim[, set_aside := outcome == "line_l1_frac" & as.character(region) ==
         REGION_LABELS[[SET_ASIDE]]]
prim[, stars := fifelse(set_aside, "", sig_stars(q))]

pA <- ggplot(prim, aes(estimate, label, colour = region, alpha = !set_aside)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    errorbar_h(aes(xmin = lo, xmax = hi),
               position = position_dodge(width = 0.66), linewidth = 0.42) +
    geom_point(size = 1.6, position = position_dodge(width = 0.66)) +
    geom_text(aes(label = stars), position = position_dodge(width = 0.66),
              hjust = -0.35, vjust = 0.75, size = 2.8, show.legend = FALSE) +
    ## The set-aside mark must ride the SAME dodge as its point, or it lands
    ## on a different region's row.
    geom_point(data = prim[set_aside == TRUE], shape = 4, size = 2.1,
               colour = PAL_CHARCOAL, show.legend = FALSE,
               position = position_dodge(width = 0.66)) +
    scale_colour_manual(values = REGION_COLORS, name = NULL) +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1), guide = "none") +
    labs(x = "Association with local SNP contribution rank (per SD)", y = NULL,
         caption = paste0(
             "\u2715 set aside: technically confounded, excluded from the ",
             "claim denominator \u2014 not a null result.\n", SIG_KEY)) +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

## ---------------------------------- b. complementary contrasts and controls
OTHER <- c("accessible_frac", "h3k27ac_frac", "h3k27me3_frac", "bivalent_frac")
oth <- prep(assoc[analysis_set == "primary" & predictor == PREDICTOR &
                  outcome %in% OTHER])
## Short, because these strips are rendered horizontally beside a facet that
## can be one row tall. The full role names are in the panel's source data and
## the caption.
ROLE_LABELS <- c(complementary_contrast = "Contrast",
                 specificity_control    = "Control",
                 descriptive            = "Descriptive")
oth[, role := factor(ROLE_LABELS[outcome_role], levels = unname(ROLE_LABELS))]

pB <- ggplot(oth, aes(estimate, label, colour = region)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    errorbar_h(aes(xmin = lo, xmax = hi),
               position = position_dodge(width = 0.66), linewidth = 0.42) +
    geom_point(size = 1.6, position = position_dodge(width = 0.66)) +
    ## Strips switched to the left and left HORIZONTAL. Rotated right-hand
    ## strips were clipped: a one-row facet has less height than "Specificity
    ## control" needs when set vertically.
    facet_grid(role ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    labs(x = "Association with local SNP contribution rank (per SD)", y = NULL,
         caption = paste0(
             "Outside the BH family of panel a; these do not share its error ",
             "rate.\nContrast = complementary contrast; Control = specificity ",
             "control.")) +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold", size = 8),
          plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"),
          plot.margin = margin(5, 10, 5, 8))

## ------------------------------------------- c. the five locked analysis sets
sens <- prep(assoc[predictor == PREDICTOR & outcome %in% BH_FAMILY])
sens[, set := factor(SET_LABELS[analysis_set], levels = unname(SET_LABELS))]
sens <- sens[!is.na(set)]
sens[, set_aside := outcome == "line_l1_frac" & as.character(region) ==
         REGION_LABELS[[SET_ASIDE]]]

pC <- ggplot(sens, aes(estimate, set, colour = region, alpha = !set_aside)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    errorbar_h(aes(xmin = lo, xmax = hi),
               position = position_dodge(width = 0.66), linewidth = 0.36) +
    geom_point(size = 1.25, position = position_dodge(width = 0.66)) +
    facet_wrap(~ flab, nrow = 1) +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1), guide = "none") +
    scale_y_discrete(limits = rev) +
    labs(x = "Association with local SNP contribution rank (per SD)", y = NULL) +
    BASE_THEME + NO_TITLES + GRID_Y

## --------------------------------- d. the continuous gradient behind panel a
##
## Deciles of the rank, not a threshold: AGENTS.md 7.2 keeps local genetic
## control continuous, and any grouping is secondary and prespecified. These
## are descriptive proportions of the same quantity panel a models.
feat <- rbindlist(lapply(regions, function(r) {
    d <- fread(file.path(rra_dir(r), "vmr-features.tsv"))
    d[, region := r][]
}), fill = TRUE)
feat[, region := as_region(region)]
feat[, decile := cut(local_snp_contribution_score, breaks = seq(0, 1, 0.1),
                     labels = 1:10, include.lowest = TRUE)]

grad <- melt(feat[!is.na(decile)],
             id.vars = c("region", "decile"),
             measure.vars = BH_FAMILY,
             variable.name = "outcome", value.name = "frac")
grad <- grad[is.finite(frac), .(mean_frac = mean(frac), n = .N),
             by = .(region, decile, outcome)]
grad[, flab := factor(OUTCOME_LABELS[as.character(outcome)],
                      levels = unname(OUTCOME_LABELS[BH_FAMILY]))]

pD <- ggplot(grad, aes(as.integer(decile), mean_frac, colour = region)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.1) +
    facet_wrap(~ flab, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_x_continuous(breaks = c(1, 5, 10)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = "Decile of local SNP contribution rank",
         y = "VMR overlap") +
    BASE_THEME + NO_TITLES

## ------------------------------------------------------------------ assemble
arm <- if (cohort == "AA") "" else paste0("_", cohort)
STEM <- paste0("figure3_repeat_repressive_architecture", arm)

figure <- (pA / pB / pC / pD) +
    plot_layout(heights = c(0.95, 1.15, 1.0, 0.85)) +
    fig_tags() & TAG_THEME

save_figure(figure, STEM, width = FIG_WIDTH_FULL, height = 8.9,
            fig_dir = fig_dir)

## ---------------------------------------------------------- source data
runs_used <- unname(RRA)
TBL <- "results/association-results-all-regions.tsv"
sd <- function(dt, nm, tbl, filt) {
    write_source_data(dt, paste0(STEM, "_", nm), runs_used, tbl, SCRIPT,
                      filt, data_dir)
}
KEEP <- c("outcome", "outcome_role", "region", "analysis_set", "predictor",
          "n", "n_fitted", "estimate", "se", "z", "p", "q")

sd(prim[, c(KEEP, "set_aside"), with = FALSE], "panel_a", TBL,
   paste0("analysis_set == 'primary'; predictor == '", PREDICTOR,
          "'; BH family (h3k9me3_frac, line_l1_frac, quiescent_frac). ",
          "set_aside marks caudate LINE/L1, which is technically confounded ",
          "and excluded from the claim denominator, NOT a null result."))
sd(oth[, KEEP, with = FALSE], "panel_b", TBL,
   paste0("analysis_set == 'primary'; predictor == '", PREDICTOR,
          "'; complementary contrasts, specificity control and descriptive ",
          "outcomes. These are outside the BH family and do not share its ",
          "error rate."))
sd(sens[, c(KEEP, "set_aside"), with = FALSE], "panel_c", TBL,
   paste0("predictor == '", PREDICTOR, "'; BH family across all five locked ",
          "analysis sets"))
sd(grad, "panel_d", "results/vmr-features.tsv",
   "all VMRs with a finite overlap fraction; deciles of local_snp_contribution_score; descriptive proportions, no enrichment test")
sd(claims, "claims", "results/interpretation-claims.tsv",
   "the accepted permitted-claim table, on the caudate gate host")

message("[done] Figure 3 written to ", fig_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
