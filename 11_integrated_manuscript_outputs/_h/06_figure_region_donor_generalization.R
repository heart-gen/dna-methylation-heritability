#### 11 / Region and donor-group generalization (Module 08) ####
##
## AGENTS.md 11, Results section 5: "shared and regionally heterogeneous
## architecture across donor groups and brain regions". Section 11 leaves
## cross-region and donor-group detail as main panels or supplements depending
## on space, so this builder names the figure by content, not by number.
##
## Every panel is one Module 08 tier, and the tier travels with it (AGENTS.md
## 7.7). Panel order follows tier rank:
##
##   A  tier 1, cross-region replication -- the direction and significance of
##      each prespecified claim-family test in each region, and the specificity
##      control that should run the other way. Direction only: no cross-region
##      comparison of levels (config agreement_basis = rank_and_direction).
##   B  tier 2, identified difference -- DLPFC minus hippocampus, both within
##      sequencing batches 1-2, over all testable pairs. A QQ plot, because the
##      question is whether regional differences exceed chance, not which test
##      is largest.
##   C  tier 3, mechanistic sensitivity -- caudate held-out R2 at its full donor
##      count and at DLPFC's n = 118, beside DLPFC. The wording is fixed by the
##      PI: "caudate magnitude attenuates after n-matching", never "inflated".
##      Tier 3 licenses whether donor count explains the excess and nothing
##      more; caudate stays batch-confounded (AGENTS.md 8.1).
##   D  donor-group axis -- AA-vs-EA ordering agreement on the pooled-discovery
##      VMR set, against the reliability ceiling it cannot exceed. Concordance
##      only: no ancestry effect, no pooled rank, no level comparison.
##
## Tier 4 (caudate vs other regions) supports nothing and has no panel; its
## estimates are the caudate column of A and the full table is Module 08's
## descriptive-confounded-regions.tsv.
##
## Usage:
##   Rscript 06_figure_region_donor_generalization.R --cohort AA --run-id ID
##   Rscript 06_figure_region_donor_generalization.R --cohort AA --run-id ID \
##       --out-dir /path/to/draft      # review draft outside the run tree

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(V2_ROOT, "11_integrated_manuscript_outputs", "_h",
                 "00_figure_theme.R"))

suppressPackageStartupMessages({
    library(patchwork)
    library(scales)
})

opts <- parse_v2_args(require = c("cohort", "run_id"))
cohort <- opts$cohort

module_root <- file.path(V2_ROOT, "11_integrated_manuscript_outputs")
out_root <- if (!is.null(opts$out_dir)) opts$out_dir else
    file.path(module_root, "_m", "runs", opts$run_id)
fig_dir  <- file.path(out_root, "figures")
data_dir <- file.path(out_root, "source_data")

## Module 08 spans all three regions; its acceptance row is keyed on the
## literal region `crossregion`.
RDG <- require_accepted_upstream("08_region_donor_generalization", cohort,
                                 "crossregion")$run_id
res_dir <- file.path(V2_ROOT, "08_region_donor_generalization", "_m", "runs",
                     RDG, "results")
rd <- function(f) {
    p <- file.path(res_dir, f)
    if (!file.exists(p)) stop("Module 08 output missing: ", p)
    as.data.table(fread(p))
}

decision <- rd("region-donor-generalization-decision.tsv")
stopifnot(decision$decision[[1]] == "PASS_REGION_DONOR_GENERALIZATION_QC")

tests       <- rd("cross-region-tests.tsv")
replication <- rd("cross-region-replication.tsv")
identified  <- rd("identified-difference.tsv")
ds_reps     <- rd("caudate-downsampling-replicates.tsv")
ds_summary  <- rd("caudate-downsampling-summary.tsv")
donor_group <- rd("donor-group-concordance.tsv")

## The interpretation constraints are guarded here, not trusted: a flag
## flipping upstream must stop the figure rather than quietly change what a
## panel is allowed to say.
stopifnot(
    all(donor_group$donor_group_inference == "concordance_only"),
    !any(as.logical(donor_group$ancestry_effect_claim_allowed)),
    !any(as.logical(donor_group$pooled_rank_emitted)),
    isTRUE(as.logical(ds_summary$caudate_remains_batch_confounded)),
    !isTRUE(as.logical(ds_summary$residual_excess_may_be_called_biological)),
    all(as.numeric(donor_group$spearman_score) <=
            as.numeric(donor_group$reliability_ceiling) + 1e-6))

rdg_cfg <- load_config("region_donor_generalization")
GROUP_LABELS <- unlist(config_get(rdg_cfg, "interpretation.group_labels"))

SCRIPT <- "11_integrated_manuscript_outputs/_h/06_figure_region_donor_generalization.R"
FIG <- paste0("figure_region_donor_generalization",
              if (cohort == "AA") "" else paste0("_", cohort))
FIG_S <- paste0("figureS_", sub("^figure_", "", FIG), "_sensitivity")

## ------------------------------------------------------------ test labels

OUTCOME_LABELS <- c(
    line_l1_frac   = "LINE/L1",
    h3k9me3_frac   = "H3K9me3",
    quiescent_frac = "Quiescent",
    h3k27me3_frac  = "H3K27me3",
    quasibinomial_continuous_local_snp_contribution_score = "CpG meQTL burden",
    expression_nearest_gene = "Nearest gene",
    expression_abc = "ABC gene",
    psi = "Splicing")
PREDICTOR_LABELS <- c(
    local_snp_contribution_score_z = "rank",
    r2_pred_oof_z         = "held-out R²",
    local_genetic_control = "rank",
    any_meqtl_support     = "any meQTL",
    meqtl_proportion      = "meQTL fraction")
FAMILY_LABELS <- c(repeat_architecture = "Repeat",
                   meqtl_burden        = "meQTL",
                   expression_coupling = "Coupling",
                   control             = "Control")

label_test <- function(outcome, predictor) {
    paste0(OUTCOME_LABELS[outcome], " · ", PREDICTOR_LABELS[predictor])
}
plotted <- replication[in_claim_family == TRUE | is_negative_control == TRUE]
for (nm in c("outcome", "predictor")) {
    lab <- if (nm == "outcome") OUTCOME_LABELS else PREDICTOR_LABELS
    unknown <- setdiff(unique(plotted[[nm]]), names(lab))
    if (length(unknown)) stop("Unlabelled ", nm, ": ", paste(unknown, collapse = ", "))
}

## Fill is direction x nominal significance, which every source module
## defines. FDR is marked separately (*), and only where the source module
## declares an FDR family: Module 04 corrects its primary prespecified tests
## only, so most repeat-architecture rows carry q = NA by design, and shading
## them "not FDR-significant" would misreport them.
STATE_FILLS <- c("Positive, P < 0.05" = "#9C4A2C",
                 "Positive, n.s."     = "#EAD2C2",
                 "Negative, n.s."     = "#D3E1EA",
                 "Negative, P < 0.05" = "#2F5F7C")
STATE_TEXT  <- c("Positive, P < 0.05" = "white", "Positive, n.s." = PAL_CHARCOAL,
                 "Negative, n.s." = PAL_CHARCOAL, "Negative, P < 0.05" = "white")

encode_state <- function(d) {
    d[, z := estimate / se]
    d[, state := factor(paste0(fifelse(estimate > 0, "Positive", "Negative"),
                               fifelse(p < 0.05, ", P < 0.05", ", n.s.")),
                        levels = names(STATE_FILLS))]
    d[, fdr_mark := fifelse(is.finite(q) & q < 0.05, "*", "")]
    d[, region := as_region(region)]
    d
}

MARK_KEY <- "* FDR < 0.05 within the source module's family   ● replicated   ○ not"

tile_theme <- BASE_THEME + NO_TITLES +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5),
          axis.text.y = element_text(size = 7.5),
          legend.position = "bottom", legend.key.size = grid::unit(0.28, "cm"),
          legend.text = element_text(size = 7),
          legend.margin = margin(-2, 0, 0, 0),
          plot.caption = element_text(size = 6.5, hjust = 0.5,
                                      colour = PAL_CHARCOAL),
          plot.caption.position = "plot")

## ------------------------------------ A. tier 1: direction in every region
##
## The claim family and the specificity controls, primary analysis set only.
## The sensitivity sets are the supplement.
rep_primary <- plotted[analysis_set == "primary"]
rep_primary[, family := fifelse(is_negative_control == TRUE, "control", analysis)]
rep_primary[, family := factor(FAMILY_LABELS[family], levels = FAMILY_LABELS)]
rep_primary[, label := label_test(outcome, predictor)]

cells <- merge(tests[analysis_set == "primary",
                     .(test_id, region, estimate, se, p, q)],
               rep_primary[, .(test_id, family, label, replicated_strict,
                               is_negative_control)],
               by = "test_id")
stopifnot(cells[, uniqueN(region), by = test_id][, all(V1 == 3L)])
cells <- encode_state(cells)

## Within each family, strongest replications first.
row_order <- cells[, .(n_sig = sum(p < 0.05), zbar = mean(z)), by = .(family, label)]
setorder(row_order, family, -n_sig, -zbar)
cells[, label := factor(label, levels = rev(unique(row_order$label)))]

rep_marks <- unique(cells[, .(family, label, replicated_strict)])
rep_marks[, `:=`(region = "Replicated",
                 mark = fifelse(replicated_strict, "●", "○"))]

pA <- ggplot(cells, aes(region, label)) +
    geom_tile(aes(fill = state), colour = "white", linewidth = 0.7) +
    geom_text(aes(label = fdr_mark, colour = state), size = 3.4, vjust = 0.78) +
    geom_text(data = rep_marks, mapping = aes(region, label, label = mark),
              size = 2.4, colour = PAL_CHARCOAL, inherit.aes = FALSE) +
    facet_grid(family ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = STATE_FILLS, name = NULL) +
    scale_colour_manual(values = STATE_TEXT, guide = "none") +
    scale_x_discrete(limits = c(REGION_ORDER, "Replicated")) +
    labs(x = NULL, y = NULL, caption = MARK_KEY) +
    guides(fill = guide_legend(ncol = 1)) +
    tile_theme +
    theme(strip.text.y = element_text(angle = -90, size = 7.5),
          panel.spacing.y = grid::unit(2.5, "pt"))

## ------------------------------- B. tier 2: DLPFC minus hippocampus, QQ plot
qq <- identified[testable == TRUE & is.finite(delta_p)]
setorder(qq, delta_p)
n_qq <- nrow(qq)
k <- seq_len(n_qq)
qq[, `:=`(expected = -log10(ppoints(n_qq)), observed = -log10(delta_p),
          lo = -log10(qbeta(0.975, k, n_qq - k + 1)),
          hi = -log10(qbeta(0.025, k, n_qq - k + 1)))]
qq[, claimed := as.logical(difference_claimed)]
qq[, label := fifelse(claimed, label_test(outcome, predictor), NA_character_)]
n_claimed <- sum(qq$claimed)

pB <- ggplot(qq, aes(expected, observed)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey90") +
    geom_abline(slope = 1, intercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    geom_point(data = qq[claimed == FALSE], colour = PAL_CHARCOAL, size = 1,
               alpha = 0.7) +
    geom_point(data = qq[claimed == TRUE], colour = PAL_RUST, size = 1.8) +
    geom_text(data = qq[claimed == TRUE], aes(label = label), hjust = 1.08,
              size = 2.3, colour = PAL_RUST) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.4, size = 2.3,
             colour = PAL_CHARCOAL, lineheight = 0.9,
             label = sprintf("DLPFC vs hippocampus\n%d tests, %d claimed",
                             n_qq, n_claimed)) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.04))) +
    labs(x = expression(Expected~-log[10]~italic(P)),
         y = expression(Observed~-log[10]~italic(P))) +
    BASE_THEME + NO_TITLES

## ---------------------- C. tier 3: caudate at its own n and at DLPFC's n
n_full <- as.integer(ds_summary$n_donors_full)
n_tgt  <- as.integer(ds_summary$target_n)
lvl <- c(sprintf("Caudate\nn = %d", n_full),
         sprintf("Caudate\nn = %d", n_tgt),
         sprintf("DLPFC\nn = %d", n_tgt))
ds_pts <- rbind(
    data.table(arm = lvl[1], region = "caudate", replicate = "full",
               mean_r2 = ds_summary$mean_r2_full),
    data.table(arm = lvl[2], region = "caudate",
               replicate = ds_reps$replicate_cell, mean_r2 = ds_reps$mean_r2_subset),
    data.table(arm = lvl[3], region = "dlpfc", replicate = "reference",
               mean_r2 = ds_summary$dlpfc_reference_mean_r2))
ds_pts[, arm := factor(arm, levels = lvl)]
ds_pts[, region := as_region(region)]
att_lab <- sprintf("%.1f%% attenuation\n95%% CI %.1f–%.1f%%",
                   100 * ds_summary$relative_attenuation_mean,
                   100 * ds_summary$relative_attenuation_ci_lower,
                   100 * ds_summary$relative_attenuation_ci_upper)
y_top <- ds_summary$mean_r2_full
y_sub <- mean(ds_reps$mean_r2_subset)

pC <- ggplot(ds_pts, aes(arm, mean_r2, colour = region)) +
    scale_x_discrete() +
    annotate("segment", x = 1.12, xend = 1.88, y = y_top, yend = y_sub + 0.001,
             colour = PAL_CHARCOAL, linewidth = 0.35,
             arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed")) +
    annotate("text", x = 1.62, y = (y_top + y_sub) / 2 + 0.005, label = att_lab,
             size = 2.3, hjust = 0, vjust = 0.5, colour = PAL_CHARCOAL,
             lineheight = 0.9) +
    geom_point(size = 2, position = position_jitter(width = 0.07, height = 0,
                                                    seed = 1)) +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0.165, 0.212), breaks = seq(0.17, 0.21, 0.01)) +
    labs(x = NULL, y = expression("Mean held-out"~R^2)) +
    BASE_THEME + NO_TITLES

## ----------------- D. donor groups: ordering agreement vs reliability ceiling
dg <- copy(donor_group)
dg[, region_lab := sprintf("%s\n%d / %d donors", REGION_LABELS[region],
                           n_donors_a, n_donors_b)]
dg[, region_lab := factor(region_lab,
                          levels = rev(region_lab[order(match(region, names(REGION_LABELS)))]))]
dg[, region_col := as_region(region)]

pD <- ggplot(dg, aes(y = region_lab)) +
    geom_segment(aes(x = 0.6, xend = reliability_ceiling, yend = region_lab),
                 colour = "grey88", linewidth = 2.4, lineend = "butt") +
    geom_point(aes(x = reliability_ceiling), shape = 124, size = 4.5,
               colour = PAL_CHARCOAL) +
    geom_point(aes(x = spearman_score, colour = region_col), size = 2.4) +
    geom_text(aes(x = spearman_score,
                  label = sprintf("%.0f%% of ceiling", 100 * fraction_of_ceiling)),
              nudge_y = 0.34, size = 2.2, colour = PAL_CHARCOAL) +
    annotate("text", x = max(dg$reliability_ceiling), y = Inf, vjust = 1.1,
             hjust = 0.5, size = 2.2, colour = PAL_CHARCOAL,
             label = "reliability\nceiling", lineheight = 0.85) +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(expand = expansion(add = c(0.5, 1))) +
    labs(x = sprintf("Spearman ρ, %s vs\n%s donors",
                     GROUP_LABELS[["AA"]], GROUP_LABELS[["EA"]]),
         y = NULL) +
    BASE_THEME + NO_TITLES +
    theme(axis.text.y = element_text(size = 7.5, lineheight = 0.9),
          axis.title.x = element_text(size = 8, lineheight = 0.95))

## ----------------------------------------------------------------- layout
right <- (pB / pC / pD) + plot_layout(heights = c(1, 0.9, 1))
figure <- (pA | right) +
    plot_layout(widths = c(0.78, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))

save_figure(figure, FIG, width = FIG_WIDTH_FULL, height = 7.4, fig_dir = fig_dir)

## ------------------------------------------------ supplement: sensitivity
##
## The claim family's repeat and repressive tests (and the H3K27me3 control)
## under every Module 04 sensitivity set: does the direction hold when
## segmental duplications, low mappability, or cell composition are handled
## differently. Same encoding as A.
SET_LABELS <- c(primary = "Primary", exclude_segdups = "Excl. segdups",
                high_mappability = "High\nmappability",
                adjust_cell_composition = "Cell-type\nadjusted",
                low_cell_composition = "Low cell-type\nvariance")
sens <- merge(tests[, .(test_id, region, estimate, se, p, q)],
              replication[analysis == "repeat_architecture" &
                              (in_claim_family == TRUE | is_negative_control == TRUE),
                          .(test_id, analysis_set, outcome, predictor,
                            is_negative_control)],
              by = "test_id")
sens <- encode_state(sens)
sens[, label := label_test(outcome, predictor)]
sens[, label := factor(label, levels = rev(unique(
    label[order(is_negative_control, outcome, predictor)])))]
sens[, set := factor(SET_LABELS[analysis_set], levels = SET_LABELS)]

pS <- ggplot(sens, aes(region, label)) +
    geom_tile(aes(fill = state), colour = "white", linewidth = 0.6) +
    geom_text(aes(label = fdr_mark, colour = state), size = 3.2, vjust = 0.78) +
    facet_grid(. ~ set) +
    scale_fill_manual(values = STATE_FILLS, name = NULL) +
    scale_colour_manual(values = STATE_TEXT, guide = "none") +
    labs(x = NULL, y = NULL,
         caption = "* FDR < 0.05 within the source module's family") +
    tile_theme +
    theme(strip.text.x = element_text(size = 7.5, face = "bold"))

save_figure(pS, FIG_S, width = FIG_WIDTH_FULL, height = 3.4, fig_dir = fig_dir)

## ------------------------------------------------------------ source data
src <- list(
    A = list(cells[, .(family, label, region, estimate, se, z, p, q, state,
                       fdr_mark, replicated_strict, is_negative_control, test_id)],
             "cross-region-tests.tsv + cross-region-replication.tsv",
             "analysis_set == primary & (in_claim_family | is_negative_control)"),
    B = list(qq[, .(test_id, analysis, outcome, predictor, delta, delta_se,
                    delta_p, delta_q, difference_claimed, expected, observed)],
             "identified-difference.tsv", "testable == TRUE"),
    C = list(ds_pts[, .(arm, region, replicate, mean_r2)],
             "caudate-downsampling-summary.tsv + caudate-downsampling-replicates.tsv",
             "all rows; DLPFC reference on DLPFC loci (no cross-region locus intersection)"),
    D = list(dg[, .(region, cell_a, cell_b, n_donors_a, n_donors_b,
                    n_loci_comparable, spearman_score, reliability_ceiling,
                    fraction_of_ceiling, donor_group_inference)],
             "donor-group-concordance.tsv", "all rows"))
for (nm in names(src)) {
    write_source_data(src[[nm]][[1]], paste0(FIG, "_", nm), RDG,
                      src[[nm]][[2]], SCRIPT, src[[nm]][[3]], data_dir)
}
write_source_data(sens[, .(analysis_set, label, region, estimate, se, z, p, q,
                           state, fdr_mark, is_negative_control, test_id)],
                  FIG_S, RDG, "cross-region-tests.tsv + cross-region-replication.tsv",
                  SCRIPT,
                  "analysis == repeat_architecture & (in_claim_family | is_negative_control)",
                  data_dir)

message("[done] ", FIG, " from ", RDG, " written to ", fig_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
