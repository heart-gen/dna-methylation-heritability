#### 11 / Figure 3: repeat-rich and repressive chromatin architecture ####
##
## AGENTS.md 11 assigns Figure 3 "LINE/L1, H3K9me3, quiescent chromatin,
## mappability, and cell sensitivities". Module 04 is the owning analysis
## (AGENTS.md 7.4); its accepted run IDs are resolved at build time.
##
## Panels, in rendered order
##   a  the BH family: H3K9me3, LINE/L1, quiescent, three regions
##   b  complementary contrasts (accessible, H3K27ac, the BrainScope ATAC
##      union) and the H3K27me3 specificity control -- the compartments that
##      run the OTHER way
##   c  the locked analysis sets, so a reader sees the sensitivities rather
##      than being told they passed
##   d  the continuous gradient behind panel a, across score deciles
##
## WHAT THIS FIGURE MAY AND MAY NOT SAY
## The permitted claim per outcome is read from the run's
## interpretation-claims.tsv at build time and asserted below; this header
## does not restate the counts, because a header that did went stale when
## Module 04 was reaccepted. Three properties the guards enforce:
##   * H3K9me3 is below the shared gate and is never drawn as shared;
##   * quiescent chromatin is shared across all regions it is required in;
##   * LINE/L1 in caudate is SET ASIDE as technically confounded -- not a
##     null, not a failure to replicate, outside the claim denominator -- and
##     renders faded with an explicit mark in panels a, c and d;
##   * a region that fails its locked sensitivities is drawn open and
##     unstarred in panels a and d, whatever its primary q.
## The seven BrainScope per-cell-type ATAC tracks are a breakdown, not a
## cell-type identification (AGENTS.md 2.3): they are carried in source data,
## outside every FDR family, and not rendered.
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
    atac_union_frac = "ATAC (union)",
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
    adjust_cell_composition = "MuSiC adjusted",
    adjust_cell_composition_scmd = "scMD adjusted",
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

## ------------------------------------------- per-region gate status, a and d
##
## The claim table carries one logical column per region: did this outcome
## survive every locked sensitivity there. A primary-model q below 0.05 is not
## that, and a region that fails its sensitivities must not be drawn as if it
## had passed -- caudate H3K9me3 has primary q 0.042 and fails high
## mappability and both composition arms. Three states, one encoding in every
## panel that draws the BH family by region:
##   supported   filled point, significance stars;
##   below_gate  open point, no stars -- estimated, not counted;
##   set_aside   faded, with a cross -- outside the claim denominator.
stopifnot(all(regions %in% names(claims)))
region_status <- function(outcome, region_key) {
    sup <- mapply(function(o, r) isTRUE(as.logical(claims[outcome == o][[r]])),
                  outcome, region_key, USE.NAMES = FALSE)
    fifelse(outcome == "line_l1_frac" & region_key %in% SET_ASIDE, "set_aside",
            fifelse(sup, "supported", "below_gate"))
}
STATUS_SHAPES <- c(supported = 16, below_gate = 21, set_aside = 16)
STATUS_ALPHA  <- c(supported = 1, below_gate = 1, set_aside = 0.35)
STATUS_LINES  <- c(supported = "solid", below_gate = "22", set_aside = "22")

## -------------------------------------------- a. the BH family, three regions
prim <- assoc[analysis_set == "primary" & predictor == PREDICTOR &
              outcome %in% BH_FAMILY]
prim[, status := region_status(outcome, region)]
prim <- prep(prim)
prim[, set_aside := status == "set_aside"]
prim[, stars := fifelse(status == "supported", sig_stars(q), "")]
prim[, mark := fifelse(set_aside, "\u2715", "")]

## The figure may not star anything the claim table does not support.
stopifnot(prim[stars != "", all(status == "supported")])
stopifnot(identical(prim[set_aside == TRUE, as.character(region)],
                    unname(REGION_LABELS[SET_ASIDE])))

## Every dodged layer is grouped by region alone. The marks are text over the
## FULL table, blank where they do not apply, so they ride the same dodge as
## their point; a layer drawn from a one-row subset is dodged on its own and
## lands on the centre row -- which once put the caudate cross on DLPFC.
DODGE <- position_dodge(width = 0.66)
pA <- ggplot(prim, aes(estimate, label, colour = region, group = region,
                       alpha = status)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    errorbar_h(aes(xmin = lo, xmax = hi), position = DODGE, linewidth = 0.42) +
    geom_point(aes(shape = status), size = 1.6, fill = "white", stroke = 0.6,
               position = DODGE) +
    geom_text(aes(label = stars), position = DODGE,
              hjust = -0.35, vjust = 0.75, size = 2.8, show.legend = FALSE) +
    geom_text(aes(label = mark), position = DODGE, colour = PAL_CHARCOAL,
              alpha = 1, size = 3.0, show.legend = FALSE) +
    scale_colour_manual(values = REGION_COLORS, name = NULL) +
    scale_shape_manual(values = STATUS_SHAPES, guide = "none") +
    scale_alpha_manual(values = STATUS_ALPHA, guide = "none") +
    labs(x = "Association with local SNP contribution rank (per SD)", y = NULL,
         caption = paste0(
             "\u25cb fails a locked sensitivity in this region: estimated, not ",
             "counted toward the claim, not starred.\n",
             "\u2715 set aside: technically confounded, excluded from the ",
             "claim denominator \u2014 not a null result.\n",
             "Stars: primary-model q, shown only where the region passes ",
             "its gate.\n", SIG_KEY)) +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

## ---------------------------------- b. complementary contrasts and controls
OTHER <- c("accessible_frac", "atac_union_frac", "h3k27ac_frac", "h3k27me3_frac",
           "bivalent_frac")
oth <- prep(assoc[analysis_set == "primary" & predictor == PREDICTOR &
                  outcome %in% OTHER])
## Short, because these strips are rendered horizontally beside a facet that
## can be one row tall. The full role names are in the panel's source data and
## the caption.
ROLE_LABELS <- c(complementary_contrast = "Contrast",
                 complementary_contrast_independent_assay = "Contrast",
                 specificity_control    = "Control",
                 descriptive            = "Descriptive")
oth[, role := factor(ROLE_LABELS[outcome_role], levels = unique(unname(ROLE_LABELS)))]
if (anyNA(oth$role)) stop("Unlabelled Module 04 outcome_role in panel b: ",
                          paste(unique(oth[is.na(role), outcome_role]), collapse = ", "))

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
## An arm Module 04 did not fit in a region (scMD: its integration gate passes
## in caudate only) is not a null and is not drawn; it is named in the caption,
## read from `arm_fitted`, so its absence cannot be mistaken for a result.
sens[, fitted := !(arm_fitted %in% c(FALSE, "FALSE")) & n_fitted > 0]
nf <- unique(sens[fitted == FALSE, .(set, region)])
nf_cap <- if (nrow(nf) == 0) NULL else paste0(
    "Not fitted, so not shown: ",
    paste(vapply(split(as.character(nf$region), as.character(nf$set)),
                 paste, character(1), collapse = ", "),
          names(split(as.character(nf$region), as.character(nf$set))),
          sep = " \u2014 ", collapse = "; "), ".")
sens <- sens[fitted == TRUE]
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
    labs(x = "Association with local SNP contribution rank (per SD)", y = NULL,
         caption = nf_cap) +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

## --------------------------------- d. the continuous gradient behind panel a
##
## Deciles of the rank, not a threshold: AGENTS.md 7.2 keeps local genetic
## control continuous, and any grouping is secondary and prespecified. These
## are descriptive proportions of the same quantity panel a models.
feat <- rbindlist(lapply(regions, function(r) {
    d <- fread(file.path(rra_dir(r), "vmr-features.tsv"))
    d[, region := r][]
}), fill = TRUE)
feat[, region_key := region]
feat[, region := as_region(region)]
feat[, decile := cut(local_snp_contribution_score, breaks = seq(0, 1, 0.1),
                     labels = 1:10, include.lowest = TRUE)]

grad <- melt(feat[!is.na(decile)],
             id.vars = c("region", "region_key", "decile"),
             measure.vars = BH_FAMILY,
             variable.name = "outcome", value.name = "frac")
grad <- grad[is.finite(frac), .(mean_frac = mean(frac), n = .N),
             by = .(region, region_key, decile, outcome)]
grad[, outcome := as.character(outcome)]
grad[, flab := factor(OUTCOME_LABELS[outcome],
                      levels = unname(OUTCOME_LABELS[BH_FAMILY]))]
## Panel d draws the same region x outcome cells as panel a, so it carries the
## same status. Drawn at full weight, caudate LINE/L1 -- set aside, primary
## estimate +0.016 -- was the steepest line in the panel and read as the
## strongest LINE/L1 result in the figure.
grad[, status := region_status(outcome, region_key)]
grad[, region_key := NULL]

pD <- ggplot(grad, aes(as.integer(decile), mean_frac, colour = region,
                       group = region, alpha = status)) +
    geom_line(aes(linetype = status), linewidth = 0.5) +
    geom_point(aes(shape = status), size = 1.1, fill = "white", stroke = 0.5) +
    facet_wrap(~ flab, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_shape_manual(values = STATUS_SHAPES, guide = "none") +
    scale_alpha_manual(values = STATUS_ALPHA, guide = "none") +
    scale_linetype_manual(values = STATUS_LINES, guide = "none") +
    scale_x_continuous(breaks = c(1, 5, 10)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = "Decile of local SNP contribution rank",
         y = "VMR overlap",
         caption = paste0(
             "Descriptive proportions; encoding as in panel a. Dashed, open: ",
             "fails a locked sensitivity in that region.\nDashed, faded: set ",
             "aside (caudate LINE/L1). Neither is counted toward the claim.")) +
    BASE_THEME + NO_TITLES +
    theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

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

sd(prim[, c(KEEP, "status", "set_aside"), with = FALSE], "panel_a", TBL,
   paste0("analysis_set == 'primary'; predictor == '", PREDICTOR,
          "'; BH family (h3k9me3_frac, line_l1_frac, quiescent_frac). ",
          "status is read from results/interpretation-claims.tsv per region: ",
          "supported (starred), below_gate (fails a locked sensitivity in ",
          "that region; not starred, not counted), set_aside (caudate ",
          "LINE/L1, technically confounded and excluded from the claim ",
          "denominator, NOT a null result)."))
sd(oth[, KEEP, with = FALSE], "panel_b", TBL,
   paste0("analysis_set == 'primary'; predictor == '", PREDICTOR,
          "'; complementary contrasts, specificity control and descriptive ",
          "outcomes. These are outside the BH family and do not share its ",
          "error rate."))
sd(sens[, c(KEEP, "set_aside"), with = FALSE], "panel_c", TBL,
   paste0("predictor == '", PREDICTOR, "'; BH family across all five locked ",
          "analysis sets"))
sd(grad, "panel_d", "results/vmr-features.tsv",
   paste0("all VMRs with a finite overlap fraction; deciles of ",
          "local_snp_contribution_score; descriptive proportions, no ",
          "enrichment test; status carried from panel a"))
ct <- prep(assoc[analysis_set == "primary" & predictor == PREDICTOR &
                 outcome_role == "celltype_breakdown_secondary"])
sd(ct[, KEEP, with = FALSE], "celltype_breakdown", TBL,
   paste0("analysis_set == 'primary'; outcome_role == 'celltype_breakdown_secondary'; ",
          "BrainScope per-cell-type ATAC tracks, outside every FDR family and NOT ",
          "rendered: a breakdown of the union contrast, not a cell-type ",
          "identification (AGENTS.md 2.3)"))
sd(claims, "claims", "results/interpretation-claims.tsv",
   "the accepted permitted-claim table, on the caudate gate host")

message("[done] Figure 3 written to ", fig_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
