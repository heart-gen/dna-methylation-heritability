#### 11 / Supplementary figures: S-LDSC, aging, environmental exposures ####
##
## Three modules whose accepted results are NULL, NOT_SUPPORTED, or exploratory
## and supplement-only. None of them is omitted on that account: a null that
## was prespecified and passed its gate is a result, and AGENTS.md 11 requires
## denominators, exclusions and the exact metric wherever a claim is made --
## including where the claim is "no effect".
##
## Figures written here
##   figureS_partitioned_heritability   Module 06, sldsc-AA-*-20260903
##   figureS_aging_axis                 Module 09b, age-AA-*-20260919
##   figureS_environmental_axis         Module 10, env-AA-*-20260920-a
##
## Each is written independently, so one blocked upstream does not stop the
## others.
##
## Usage:
##   Rscript 11_supplementary_figures.R --cohort AA --run-id fig-all-YYYYMMDD
##   Rscript 11_supplementary_figures.R --cohort AA --run-id smoke \
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

SCRIPT <- "11_integrated_manuscript_outputs/_h/11_supplementary_figures.R"
arm <- if (cohort == "AA") "" else paste0("_", cohort)

by_region <- function(fun) rbindlist(lapply(regions, function(r) {
    d <- fun(r); if (is.null(d) || nrow(d) == 0) return(NULL); d[, region := r][]
}), fill = TRUE)

## ===================================================== S-LDSC (Module 06)
##
## The accepted result is a NULL: sldsc_supports_brain_enrichment = FALSE in
## all three cells, 0 of 8 traits FDR-significant. The figure has to make the
## null legible rather than imply a signal, so it plots the enrichment estimate
## WITH its standard error -- the SEs are large, which is the actual finding --
## and states the FDR outcome on the panel.
SLDSC <- vapply(regions, function(r)
    require_accepted_upstream("06_partitioned_heritability", cohort, r)$run_id,
    character(1))
sldsc_dir <- function(r) file.path(V2_ROOT, "06_partitioned_heritability",
                                   "_m", "runs", SLDSC[[r]], "results")

met <- by_region(function(r) fread(file.path(sldsc_dir(r), "sldsc-metrics.tsv")))
dec <- by_region(function(r) fread(file.path(sldsc_dir(r), "partitioned-h2-decision.tsv")))

## The null is the accepted reading; if a cell ever turns positive the panel's
## framing is wrong.
supports <- grep("supports", names(dec), value = TRUE)
if (length(supports) == 1 && any(dec[[supports]] %in% c(TRUE, "TRUE"))) {
    stop("A Module 06 cell now supports brain enrichment; this supplement is ",
         "written as a null. Re-read the accepted decision.")
}
met[, region := as_region(region)]
met[, `:=`(lo = enrichment - 1.96 * enrichment_se,
           hi = enrichment + 1.96 * enrichment_se)]
met[, class := factor(fifelse(trait_class == "brain", "Brain", "Non-brain control"),
                      levels = c("Brain", "Non-brain control"))]
## Order by the pooled estimate so the panel reads as one ranking.
ord <- met[, .(m = mean(enrichment)), by = trait_label][order(-m)]
met[, trait_label := factor(trait_label, levels = ord$trait_label)]
n_sig <- met[, sum(tau_q < 0.05, na.rm = TRUE)]

pS1 <- ggplot(met, aes(enrichment, trait_label, colour = region)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    errorbar_h(aes(xmin = lo, xmax = hi),
               position = position_dodge(width = 0.7), linewidth = 0.4) +
    geom_point(size = 1.4, position = position_dodge(width = 0.7)) +
    facet_grid(class ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_manual(values = REGION_COLORS, name = NULL) +
    labs(x = "S-LDSC enrichment of the local genetic-control annotation",
         y = NULL,
         caption = paste0(
             "Prespecified 8-trait family, EUR LD scores. ", n_sig,
             " of ", nrow(met), " tests reach q < 0.05: the annotation shows",
             "\nno partitioned-heritability enrichment in any region. Wide",
             " intervals are the result, not a rendering artefact.\nThe",
             " annotation is a genomic feature, not a donor-group LD claim.")) +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold", size = 8),
          plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

S1 <- paste0("figureS_partitioned_heritability", arm)
save_figure(pS1, S1, width = FIG_WIDTH_FULL, height = 4.2, fig_dir = fig_dir)
write_source_data(met[, .(region, trait_label, trait_class, annotation,
                          prop_snps, prop_h2, prop_h2_se, enrichment,
                          enrichment_se, enrichment_p, tau, tau_se, tau_z,
                          tau_p, tau_q, total_h2, total_h2_se, lambda_gc,
                          intercept)],
                  paste0(S1, "_panel_a"), unname(SLDSC),
                  "results/sldsc-metrics.tsv", SCRIPT,
                  "the frozen 8-trait family; accepted decision is a NULL (sldsc_supports_brain_enrichment = FALSE in all three cells)",
                  data_dir)

## ===================================================== Aging (Module 09b)
##
## Cross-region token is NOT_SUPPORTED. Direction is concordant in all three
## regions and two are individually significant, but the cell_composition_r2
## gating arm removes the gradient everywhere -- so the sensitivity arms are a
## PANEL, not a footnote. The reading is "the age-responsive low-control VMRs
## are the composition-sensitive ones", and the figure must not be readable as
## a clean aging result.
AGE <- vapply(regions, function(r)
    require_accepted_upstream("09b_aging_application", cohort, r)$run_id,
    character(1))
age_dir <- function(r) file.path(V2_ROOT, "09b_aging_application", "_m",
                                 "runs", AGE[[r]], "results")

per <- fread(file.path(V2_ROOT, "09b_aging_application", "_m", "combined",
                       sprintf("aging-per-region-%s.tsv", cohort)))
gat <- by_region(function(r) fread(file.path(age_dir(r), "gating-sensitivities.tsv")))
qrt <- by_region(function(r) fread(file.path(age_dir(r), "axis-quartile-summary.tsv")))

stopifnot(all(per$cross_sectional_design %in% c(TRUE, "TRUE")))
stopifnot(all(per$causal_interpretation_allowed %in% c(FALSE, "FALSE")))
per[, region := as_region(region)]
gat[, region := as_region(region)]

pA1 <- ggplot(per, aes(primary_estimate, region, colour = region)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    errorbar_h(aes(xmin = primary_ci_lower, xmax = primary_ci_upper),
               linewidth = 0.45) +
    geom_point(size = 1.8) +
    geom_text(aes(label = paste0("P = ", signif(primary_p, 2))),
              x = Inf, hjust = 1.05, vjust = -1.0, size = 2.2,
              colour = "grey35") +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_y_discrete(limits = rev, expand = expansion(add = c(0.5, 0.8))) +
    labs(x = "Primary: debiased squared age effect per SD of rank", y = NULL,
         caption = paste("Cross-sectional design: age-associated",
                         "differences,\nnever change with age.")) +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

## The arm that removes the gradient gets the same visual weight as the primary.
MEMBER_LABELS <- c(controls_only = "Controls only",
                   cell_music = "MuSiC cell PCs",
                   cell_scmd = "scMD cell PCs",
                   cell_composition_r2 = "VMR composition sensitivity")
gat[, mlab := factor(MEMBER_LABELS[member], levels = rev(unname(MEMBER_LABELS)))]
gat <- gat[!is.na(mlab)]

## A member that was never FITTED is not a null, and ggplot dropping it with
## "Removed 2 rows containing missing values" leaves a reader unable to tell
## the two apart. scMD PCs are fitted only where the DNAm integration gate
## passes, which is caudate alone (AGENTS.md 7.9), so the other two cells are
## named in the caption rather than silently vanishing.
not_fitted <- gat[!(fitted %in% c(TRUE, "TRUE"))]
nf_note <- if (nrow(not_fitted) == 0) "" else paste0(
    "\nNot fitted, so not shown: ",
    paste(sprintf("%s in %s", MEMBER_LABELS[not_fitted$member],
                  as.character(not_fitted$region)), collapse = "; "),
    " (the scMD integration gate passes in caudate only).")
gat_fit <- gat[fitted %in% c(TRUE, "TRUE")]

pA2 <- ggplot(gat_fit, aes(estimate, mlab, colour = region,
                           shape = survives)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    geom_point(size = 1.7, position = position_dodge(width = 0.66)) +
    scale_colour_manual(values = REGION_COLORS, name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4),
                       labels = c(`TRUE` = "survives", `FALSE` = "fails"),
                       name = NULL) +
    labs(x = "Gating-sensitivity estimate", y = NULL,
         caption = paste0("The VMR composition-sensitivity arm removes the ",
                          "gradient in every region, which is why\nthe ",
                          "cross-region token is NOT_SUPPORTED. Read as: the ",
                          "age-responsive low-control\nVMRs are the ",
                          "composition-sensitive ones.", nf_note)) +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

S2 <- paste0("figureS_aging_axis", arm)
save_figure((pA1 / pA2) + plot_layout(heights = c(0.8, 1)) + fig_tags() & TAG_THEME,
            S2, width = FIG_WIDTH_THREEQ, height = 6.0, fig_dir = fig_dir)
write_source_data(per, paste0(S2, "_panel_a"), unname(AGE),
                  "_m/combined/aging-per-region-{cohort}.tsv", SCRIPT,
                  "accepted per-region rows; cross-sectional design; causal interpretation not allowed",
                  data_dir)
write_source_data(gat, paste0(S2, "_panel_b"), unname(AGE),
                  "results/gating-sensitivities.tsv", SCRIPT,
                  "all gating members under strict conjunction; cell_composition_r2 fails in every region",
                  data_dir)
write_source_data(qrt, paste0(S2, "_quartiles"), unname(AGE),
                  "results/axis-quartile-summary.tsv", SCRIPT,
                  "secondary quartile contrast; the continuous model is primary (AGENTS.md 7.2)",
                  data_dir)

## =============================================== Environmental (Module 10)
##
## Exploratory, supplement-only, and gated on acceptance (AGENTS.md 6). If the
## acceptance rows are not in the Module 10 README this block is SKIPPED with a
## message rather than failing the run -- the other supplements do not depend
## on it.
env_ok <- tryCatch({
    invisible(vapply(regions, function(r)
        require_accepted_upstream("10_environmental_exploratory", cohort, r)$run_id,
        character(1)))
    TRUE
}, error = function(e) {
    message("[skip] figureS_environmental_axis: ", conditionMessage(e))
    FALSE
})

if (env_ok) {
    ecomb <- file.path(V2_ROOT, "10_environmental_exploratory", "_m", "combined")
    eax <- fread(file.path(ecomb, sprintf("environmental-axis-per-region-%s.tsv", cohort)))
    ever <- fread(file.path(ecomb, sprintf("environmental-vmr-associations-fdr-%s.tsv", cohort)))

    if (!all(eax$citable %in% c(TRUE, "TRUE"))) {
        stop("Module 10's collated tables still carry citable = FALSE. ",
             "Re-run 10/_h/06_collate_regions.R WITHOUT --allow-unaccepted-runs.")
    }
    stopifnot(all(eax$exploratory_supplement_only %in% c(TRUE, "TRUE")))
    stopifnot(all(eax$environmentally_determined_claim_allowed %in% c(FALSE, "FALSE")))

    eax[, region := as_region(region)]
    eax[, fam := paste0(exposure, " @ ", stratum)]
    ## Ratio-scale families only: a `raw` family has mean(omega) <= 0, so no
    ## ratio exists and its estimate is not on the same scale.
    rel <- eax[primary_scale != "raw"]
    rel[, fam := factor(fam, levels = unique(rel[order(primary_beta), fam]))]
    rel[, stars := sig_stars(primary_fdr)]

    ## One null family (dlpfc marital_status, the family the retired
    ## relative-scale guard was written for) has a point estimate near -2.4 and
    ## an interval wider than the rest of the figure put together. On a common
    ## scale it squeezes the three FDR-surviving families into a sliver. The
    ## window below is STATED, and anything outside it is drawn at the boundary
    ## with its real numbers printed, so nothing is hidden by the clip.
    LIM <- c(-1.3, 0.7)
    rel[, offscale := primary_ci_lower < LIM[1] | primary_ci_upper > LIM[2] |
            primary_beta < LIM[1] | primary_beta > LIM[2]]
    rel[, off_lab := sprintf("%.2f [%.2f, %.2f]", primary_beta,
                             primary_ci_lower, primary_ci_upper)]

    pE <- ggplot(rel, aes(primary_beta, fam, colour = region)) +
        geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
        errorbar_h(aes(xmin = primary_ci_lower, xmax = primary_ci_upper),
                   position = position_dodge(width = 0.72), linewidth = 0.4) +
        geom_point(size = 1.5, position = position_dodge(width = 0.72)) +
        geom_text(aes(label = stars), position = position_dodge(width = 0.72),
                  hjust = -0.4, vjust = 0.75, size = 2.5, show.legend = FALSE) +
        geom_text(data = rel[offscale == TRUE], aes(label = off_lab),
                  x = LIM[1], hjust = -0.03, vjust = -0.9, size = 2.0,
                  show.legend = FALSE) +
        coord_cartesian(xlim = LIM) +
        scale_colour_manual(values = REGION_COLORS, name = NULL) +
        labs(x = "Proportional gradient in exposure-explained variance per SD of rank",
             y = NULL,
             caption = paste0(
                 "Exploratory supplement, x axis clipped to [", LIM[1], ", ",
                 LIM[2], "]; estimates outside it are printed in full.\n",
                 "NO percentage here is formally identifiable: the family mean ",
                 "never separates from zero (max |z| = ",
                 sprintf("%.2f", max(eax$mean_omega_z, na.rm = TRUE)),
                 " < 1.96),\nso the ratio replicates but the level it is a ",
                 "ratio of does not. The donor bootstrap inflates the ",
                 "denominator ",
                 sprintf("%.1f", min(eax$bootstrap_mean_omega_inflation, na.rm = TRUE)),
                 "-",
                 sprintf("%.1f", max(eax$bootstrap_mean_omega_inflation, na.rm = TRUE)),
                 "x,\nso these p-values may be too small. A negative gradient ",
                 "is NOT evidence that exposure effects concentrate in weakly ",
                 "controlled VMRs.")) +
        BASE_THEME + NO_TITLES + GRID_Y +
        theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
              plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

    S3 <- paste0("figureS_environmental_axis", arm)
    save_figure(pE, S3, width = FIG_WIDTH_FULL, height = 4.0, fig_dir = fig_dir)
    write_source_data(eax, paste0(S3, "_panel_a"),
                      unique(eax$run_id),
                      "_m/combined/environmental-axis-per-region-{cohort}.tsv",
                      SCRIPT,
                      "all stage B families; the rendered panel shows ratio-scale families only, because a `raw` family has mean(omega) <= 0 and no ratio exists",
                      data_dir)
    write_source_data(ever, paste0(S3, "_stage_a"),
                      unique(eax$run_id),
                      "_m/combined/environmental-vmr-associations-fdr-{cohort}.tsv",
                      SCRIPT,
                      "stage A FDR-surviving VMR x exposure pairs (0 caudate / 0 dlpfc / 12 hippocampus)",
                      data_dir)
}

## ============================================ Schizophrenia (Module 09)
##
## The schizophrenia locus evidence, displaced here from Figure 5. Module 09's
## decision 2 is scz_application_retention = RETAIN_MAIN_TEXT and that is
## honoured in Figure 5, where schizophrenia appears as a marked example among
## the other 62 traits. What belongs in the supplement is the locus-level
## detail, which is specific to this trait and would otherwise crowd out the
## trait-general result the main figure exists to show.
##
## The direction is a DEPLETION: SCZ-linked VMRs sit LOWER on the axis in all
## three regions. It is never rendered as an enrichment (AGENTS.md 7.8).
SCZR <- vapply(regions, function(r)
    require_accepted_upstream("09_schizophrenia_risk_application", cohort, r)$run_id,
    character(1))
scz_dir <- function(r) file.path(V2_ROOT, "09_schizophrenia_risk_application",
                                 "_m", "runs", SCZR[[r]], "results")

axt <- by_region(function(r) fread(file.path(scz_dir(r), "architecture-axis-tests.tsv")))
lev <- by_region(function(r) fread(file.path(scz_dir(r), "locus-evidence.tsv")))

## The result is a depletion in every region. If that ever flips, this figure's
## axis labels and caption are wrong.
prim <- axt[predictor == "local_snp_contribution_score_z" &
            model == "wilcoxon_rank_sum"]
if (!all(prim$direction == "lower_in_scz_linked")) {
    stop("Module 09's axis contrast is no longer lower_in_scz_linked in every ",
         "region; this supplement is written as a depletion.")
}
if (any(grepl("enrich", prim$contrast_label, ignore.case = TRUE) &
        !grepl("DEPLETION", prim$contrast_label))) {
    stop("Module 09 emitted an enrichment framing; AGENTS.md 7.8 forbids it.")
}
prim[, region := as_region(region)]
lev[, region := as_region(region)]

pZ1 <- ggplot(prim, aes(estimate, region, colour = region)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    geom_point(size = 2.1) +
    geom_text(aes(label = paste0(sig_stars(qvalue), "  n = ",
                                 label_comma()(n_linked), " linked / ",
                                 label_comma()(n_background), " background")),
              x = Inf, hjust = 1.03, vjust = -1.0, size = 2.2,
              colour = "grey35") +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_y_discrete(limits = rev, expand = expansion(add = c(0.5, 0.9))) +
    labs(x = "Mean axis difference, SCZ-linked minus background (score SD)",
         y = NULL,
         caption = paste("Negative = SCZ-linked VMRs sit LOWER on the",
                         "local genetic-control axis. A depletion, never an",
                         "enrichment.\nThe depletion is trait-general (Fig. 5);",
                         "schizophrenia is a typical example, not a",
                         "schizophrenia-specific effect.")) +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

## Locus evidence tiers. Counts of loci carrying each independent line of
## support, so the reader sees how thin or thick each tier is.
EV <- c(has_significant_risk_variant_cpg_meqtl = "CpG meQTL support",
        has_transcriptional_coupling          = "Transcriptional coupling",
        has_external_genetic_support          = "GTEx support",
        has_claimable_colocalization          = "Claimable colocalization")
ev <- rbindlist(lapply(names(EV), function(k) {
    lev[, .(tier = EV[[k]], n_loci = sum(get(k) %in% c(TRUE, "TRUE")),
            n_total = .N), by = region]
}))
ev[, tier := factor(tier, levels = rev(unname(EV)))]

pZ2 <- ggplot(ev, aes(n_loci, tier, fill = region)) +
    geom_col(width = 0.7, position = position_dodge(width = 0.76)) +
    geom_text(aes(label = n_loci), position = position_dodge(width = 0.76),
              hjust = -0.25, size = 2.2, colour = "black") +
    scale_fill_manual(values = REGION_COLORS, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = "Schizophrenia-risk loci with the evidence type", y = NULL,
         caption = paste0("Out of ", lev[, .N, by = region][, paste(N, collapse = "/")],
                          " loci linked to a VMR (caudate/DLPFC/hippocampus).",
                          "\nGWAS loci defined from European-ancestry summary",
                          " statistics; the cohort is admixed African American.",
                          "\nCaudate carries CAUDATE_MAGNITUDE_CLAIM_NOT_SUPPORTED",
                          " and is batch-confounded.")) +
    BASE_THEME + NO_TITLES + GRID_X +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

## The prioritized loci, as an evidence grid rather than a table of ticks.
pri <- lev[prioritized %in% c(TRUE, "TRUE")]
pri_long <- rbindlist(lapply(names(EV), function(k) {
    pri[, .(region, locus = paste0(chrom, ":", index_snp), tier = EV[[k]],
            has = get(k) %in% c(TRUE, "TRUE"))]
}))
pri_long[, tier := factor(tier, levels = unname(EV))]

pZ3 <- ggplot(pri_long, aes(tier, locus, fill = has)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    facet_wrap(~ region, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c(`TRUE` = PAL_RUST, `FALSE` = "#EFEAE4"),
                      labels = c(`TRUE` = "present", `FALSE` = "absent"),
                      name = NULL) +
    labs(x = NULL, y = NULL,
         caption = paste("Up to five prioritized loci per region under the",
                         "prespecified ranked composite rule.\nAn elastic-net",
                         "or index SNP is not a causal variant, and no panel",
                         "here claims mediation.")) +
    BASE_THEME + NO_TITLES +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7),
          plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

S4 <- paste0("figureS_schizophrenia_application", arm)
save_figure((pZ1 / pZ2 / pZ3) + plot_layout(heights = c(0.8, 1.0, 1.15)) +
                fig_tags() & TAG_THEME,
            S4, width = FIG_WIDTH_FULL, height = 8.6, fig_dir = fig_dir)

write_source_data(prim, paste0(S4, "_panel_a"), unname(SCZR),
                  "results/architecture-axis-tests.tsv", SCRIPT,
                  "predictor == 'local_snp_contribution_score_z'; model == 'wilcoxon_rank_sum'; direction is lower_in_scz_linked (a DEPLETION) in all three regions; trait-general, never schizophrenia-specific",
                  data_dir)
write_source_data(ev, paste0(S4, "_panel_b"), unname(SCZR),
                  "results/locus-evidence.tsv", SCRIPT,
                  "counts of SCZ-risk loci carrying each evidence type; GWAS loci are European-ancestry; caudate is batch-confounded",
                  data_dir)
write_source_data(pri, paste0(S4, "_panel_c"), unname(SCZR),
                  "results/locus-evidence.tsv", SCRIPT,
                  "prioritized == TRUE; prespecified ranked composite rule, at most five loci per region",
                  data_dir)

message("[done] supplementary figures written to ", fig_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
