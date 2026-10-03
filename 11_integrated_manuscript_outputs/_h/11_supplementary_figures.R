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
## The accepted result is a NULL under the two-annotation model
## (sldsc-AA-*-20260925): tau of LOCAL_SNP_CONTRIBUTION_Z is estimated
## CONDITIONAL on VMR_TESTED, so it is the within-VMR gradient, and no trait
## reaches q < 0.05 in any cell.
##
## What is plotted is that tau, as a z-score. Until 2026-10-03 this panel
## plotted `enrichment`, which for a signed continuous annotation is a ratio
## over a signed sum and not a share of heritability; Module 06 marks exactly
## those rows enrichment_interpretable = FALSE. The guard below stops the build
## if an uninterpretable enrichment would be drawn.
##
## Module 06's writing rule: "no detectable enrichment at this footprint",
## never "no enrichment" -- the annotation is small and the module has no
## positive control, so a null cannot be told from a power null.
SLDSC <- vapply(regions, function(r)
    require_accepted_upstream("06_partitioned_heritability", cohort, r)$run_id,
    character(1))
sldsc_dir <- function(r) file.path(V2_ROOT, "06_partitioned_heritability",
                                   "_m", "runs", SLDSC[[r]], "results")

met <- by_region(function(r) fread(file.path(sldsc_dir(r), "sldsc-metrics.tsv")))
dec <- by_region(function(r) fread(file.path(sldsc_dir(r), "partitioned-h2-decision.tsv")))

if (!all(dec$tau_conditional_on_vmr_membership %in% c(TRUE, "TRUE"))) {
    stop("A Module 06 cell is not the two-annotation model; its tau is not the ",
         "within-VMR gradient this panel describes.")
}
supports <- grep("supports", names(dec), value = TRUE)
if (length(supports) == 1 && any(dec[[supports]] %in% c(TRUE, "TRUE"))) {
    stop("A Module 06 cell now supports brain enrichment; this supplement is ",
         "written as a null. Re-read the accepted decision.")
}
met <- met[is_primary_hypothesis %in% c(TRUE, "TRUE")]
if (any(met$enrichment_interpretable %in% c(TRUE, "TRUE"))) {
    stop("A primary-hypothesis S-LDSC row is now marked enrichment-interpretable; ",
         "re-decide what this panel should show.")
}
met[, region := as_region(region)]
met[, class := factor(fifelse(trait_class == "brain", "Brain", "Non-brain control"),
                      levels = c("Brain", "Non-brain control"))]
## The facet already says "Non-brain control"; the label suffix repeating it
## pushed the panel to a sliver.
met[, trait_label := sub("\\s*--\\s*non-brain negative control$", "", trait_label)]
ord <- met[, .(m = mean(tau_z)), by = trait_label][order(-m)]
met[, trait_label := factor(trait_label, levels = ord$trait_label)]
n_sig <- met[, sum(tau_q < 0.05, na.rm = TRUE)]
n_traits <- uniqueN(met$trait_label)

pS1 <- ggplot(met, aes(tau_z, trait_label, colour = region)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    geom_vline(xintercept = c(-1.96, 1.96), colour = PAL_NULL, linewidth = 0.3,
               linetype = 2) +
    geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
    facet_grid(class ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_manual(values = REGION_COLORS, name = NULL) +
    labs(x = "S-LDSC \u03c4 z, score | VMR membership",
         y = NULL,
         caption = paste(strwrap(paste0(
             "Prespecified ", n_traits, "-trait family, EUR LD scores, two-annotation ",
             "model. ", n_sig, " of ", nrow(met), " tests reach q < 0.05: no detectable ",
             "enrichment at this footprint. The annotation is small and there is no ",
             "positive control, so this null is not distinguishable from a power null. ",
             "Dashed lines: |z| = 1.96 (nominal)."), width = 115), collapse = "\n")) +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, face = "bold", size = 8),
          plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"),
          plot.caption.position = "plot")

S1 <- paste0("figureS_partitioned_heritability", arm)
save_figure(pS1, S1, width = FIG_WIDTH_FULL, height = 4.2, fig_dir = fig_dir)
write_source_data(met[, .(region, trait_label, trait_class, annotation_name,
                          annotation_role, is_primary_hypothesis,
                          enrichment_interpretable, tau, tau_se, tau_z, tau_p,
                          tau_q, in_fdr_family, n_tests_in_family, total_h2,
                          total_h2_se, lambda_gc, intercept, ld_reference_arm)],
                  paste0(S1, "_panel_a"), unname(SLDSC),
                  "results/sldsc-metrics.tsv", SCRIPT,
                  paste0("is_primary_hypothesis (LOCAL_SNP_CONTRIBUTION_Z, tau conditional on VMR_TESTED); ",
                         n_sig, " of ", nrow(met),
                         " q < 0.05; enrichment omitted because enrichment_interpretable = FALSE"),
                  data_dir)

## ===================================================== Aging (Module 09b)
##
## The design constraint, which does not depend on the result: a gating arm
## that removes the gradient is the whole reading of this module, so the
## sensitivity arms get a PANEL rather than a footnote, and the figure must
## not be readable as a clean aging result whatever the arms do.
##
## The VERDICT -- the cross-region token, which regions are supported, which
## arms fail -- is read out of the run's own decision tables below and never
## written into a string here. A caption that states a per-region verdict as
## prose goes stale silently the next time the verdict changes, and this one
## did: it asserted "removes the gradient in every region, which is why the
## cross-region token is NOT_SUPPORTED" while the module 09b scMD-gate
## correction was in flight. AGENTS.md 7.11 wants a panel's numbers traceable
## to a run, table, script and filter; a hard-coded claim is none of those.
AGE <- vapply(regions, function(r)
    require_accepted_upstream("09b_aging_application", cohort, r)$run_id,
    character(1))
age_dir <- function(r) file.path(V2_ROOT, "09b_aging_application", "_m",
                                 "runs", AGE[[r]], "results")

per <- fread(file.path(V2_ROOT, "09b_aging_application", "_m", "combined",
                       sprintf("aging-per-region-%s.tsv", cohort)))
gat <- by_region(function(r) fread(file.path(age_dir(r), "gating-sensitivities.tsv")))
qrt <- by_region(function(r) fread(file.path(age_dir(r), "axis-quartile-summary.tsv")))

## The stage-05 cross-region decision row, which is where the token lives.
xr <- fread(file.path(V2_ROOT, "09b_aging_application", "_m", "combined",
                      sprintf("aging-cross-region-decision-%s.tsv", cohort)))
if (nrow(xr) != 1L) {
    stop("Expected one aging cross-region decision row for ", cohort,
         ", found ", nrow(xr))
}
XR_TOKEN <- as.character(xr$aging_axis_association)
if (!nzchar(XR_TOKEN) || is.na(XR_TOKEN)) {
    stop("Aging cross-region decision carries no aging_axis_association token.")
}

stopifnot(all(per$cross_sectional_design %in% c(TRUE, "TRUE")))
stopifnot(all(per$causal_interpretation_allowed %in% c(FALSE, "FALSE")))
## The 2026-10-03 magnitude gate (T28 decision B): a region whose family mean
## is not separated from zero reports a SIGNED TEST only, never a proportional
## change. The label says which, from the run's own column, so a reader cannot
## take the point estimate of a gated region as a percentage.
per[, mag_ok := primary_relative_magnitude_reportable %in% c(TRUE, "TRUE")]
per[, p_lab := paste0("P = ", signif(primary_p, 2),
                      fifelse(mag_ok, "", " \u00b7 signed test only"))]
per[, region := as_region(region)]
gat[, region := as_region(region)]

pA1 <- ggplot(per, aes(primary_estimate, region, colour = region)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    errorbar_h(aes(xmin = primary_ci_lower, xmax = primary_ci_upper),
               linewidth = 0.45) +
    geom_point(size = 1.8) +
    geom_text(aes(label = p_lab),
              x = Inf, hjust = 1.05, vjust = -1.0, size = 2.2,
              colour = "grey35") +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_y_discrete(limits = rev, expand = expansion(add = c(0.5, 0.8))) +
    labs(x = "Primary gradient per SD of rank", y = NULL,
         caption = paste(strwrap(paste0(
             "Gradient in the debiased squared age effect, ratio scale. ",
             "Cross-sectional design: age-associated differences, never change ",
             "with age. \"Signed test only\": the magnitude gate withholds a ",
             "proportional reading (", sum(!per$mag_ok), " of ", nrow(per),
             " regions)."), width = 78), collapse = "\n")) +
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
## Module 09b's scMD-gate correction adds a `reason` column, so the reason can
## be quoted from the run instead of assumed. Before it lands there is no such
## column and the note simply says the arm was not fitted -- which is the point:
## the caption never asserts why on the run's behalf.
nf_reason <- if ("reason" %in% names(not_fitted))
    as.character(not_fitted$reason) else rep(NA_character_, nrow(not_fitted))
nf_note <- if (nrow(not_fitted) == 0) "" else paste0(
    "\nNot fitted, so not shown: ",
    paste(sprintf("%s in %s%s", MEMBER_LABELS[not_fitted$member],
                  as.character(not_fitted$region),
                  ifelse(is.na(nf_reason) | !nzchar(nf_reason), "",
                         paste0(" (", nf_reason, ")"))),
          collapse = "; "), ".")
gat_fit <- gat[fitted %in% c(TRUE, "TRUE")]

## ------------------------------------------- panel b caption, derived not typed
##
## Every quantity below comes out of the tables already loaded. "Every region"
## is a claim about a denominator, so the denominator is counted; which regions
## are supported is read from `region_supported`; the token is read from the
## stage-05 decision row. Nothing here changes if the verdict changes -- only
## what it prints.
failing <- gat_fit[!(survives %in% c(TRUE, "TRUE"))]
fail_note <- if (nrow(failing) == 0) {
    "No fitted gating arm removes the gradient."
} else {
    tab <- merge(failing[, .(n_fail = .N), by = member],
                 gat_fit[, .(n_fitted = .N), by = member], by = "member")
    paste0("Gradient removed by: ",
           paste(sprintf("%s (%d of %d regions fitted)",
                         MEMBER_LABELS[tab$member], tab$n_fail, tab$n_fitted),
                 collapse = "; "), ".")
}
supported <- per[region_supported %in% c(TRUE, "TRUE"), as.character(region)]
support_note <- if (length(supported) == 0)
    "No region's reading is supported." else
    paste0("Supported in: ", paste(sort(supported), collapse = ", "), ".")
## No interpretive sentence is attached to a failing arm. The one that used to
## be here -- "the age-responsive low-control VMRs are the arm-sensitive ones"
## -- is the composition qualifier the 2026-10-01 reacceptance retired: it was
## true only of the scMD-derived cell_composition_r2 arm, which Module 04's
## 2026-09-25 reacceptance replaced with a MuSiC-derived one.
## A derived caption has no fixed length, so it is wrapped to the panel rather
## than hand-broken. The hard-coded one could carry its own newlines because
## nobody expected it to change.
CAP_B <- paste(unlist(lapply(
    c(paste0("Cross-region token: ", XR_TOKEN, ". ", support_note),
      fail_note,
      sub("^\n", "", nf_note)),
    function(s) if (!nzchar(s)) NULL else strwrap(s, width = 78))),
    collapse = "\n")

pA2 <- ggplot(gat_fit, aes(estimate, mlab, colour = region,
                           shape = survives)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    geom_point(size = 1.7, position = position_dodge(width = 0.66)) +
    scale_colour_manual(values = REGION_COLORS, name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 4),
                       labels = c(`TRUE` = "survives", `FALSE` = "fails"),
                       name = NULL) +
    labs(x = "Gating-sensitivity estimate", y = NULL, caption = CAP_B) +
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
## `row_filter` describes the rows, not the result. The verdict the old string
## asserted -- which arm fails where -- is carried as DATA here: `fitted` and
## `survives` per member per region, plus the region reading and the stage-05
## token the caption prints. A reader checking the caption against the table
## finds the same values, because the caption is built from them.
gat_sd <- merge(gat,
                per[, .(region, region_reading, region_supported)],
                by = "region", all.x = TRUE)
gat_sd[, `:=`(cross_region_token = XR_TOKEN, caption_rendered = CAP_B)]
write_source_data(gat_sd, paste0(S2, "_panel_b"), unname(AGE),
                  "results/gating-sensitivities.tsv + _m/combined/aging-cross-region-decision-{cohort}.tsv",
                  SCRIPT,
                  "all gating members under strict conjunction; not-fitted members retained and flagged, never dropped",
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

    ## PI 2026-10-03: one family's FDR call is undetermined across bootstrap
    ## seeds and must not be drawn as resolved either way. Matched before the
    ## region labels are prettified, against the raw region token.
    eax[, fdr_call_undetermined := apply_reporting_constraint(
        read_reporting_constraints(run_dir), "fdr_call_undetermined",
        "environmental-axis-per-region", eax, "region",
        paste0(eax$exposure, "@", eax$stratum))]
    n_fam <- nrow(eax)
    n_reportable <- sum(eax$relative_magnitude_reportable %in% c(TRUE, "TRUE"))
    max_den_z <- max(as.numeric(eax$relative_magnitude_den_z), na.rm = TRUE)
    eax[, region := as_region(region)]
    eax[, fam := paste0(exposure, " @ ", stratum)]
    ## Ratio-scale families only: a `raw` family has mean(omega) <= 0, so no
    ## ratio exists and its estimate is not on the same scale.
    rel <- eax[primary_scale != "raw"]
    rel[, fam := factor(fam, levels = unique(rel[order(primary_beta), fam]))]
    rel[, stars := data.table::fifelse(fdr_call_undetermined, "?",
                                       sig_stars(primary_fdr))]

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

    ## Derived, so of no fixed length: wrapped to the panel, never hand-broken
    ## (ggplot does not wrap captions; the unwrapped draft ran off the page).
    CAP_E <- paste(unlist(lapply(c(
        paste0("Exploratory supplement; x axis clipped to [", LIM[1], ", ", LIM[2],
               "], estimates outside it printed in full. ", SIG_KEY,
               "   ? FDR call undetermined across bootstrap seeds."),
        if (n_reportable == 0L) paste0(
            "Signed gradients only: no family passes the relative-magnitude gate (0 of ",
            n_fam, "; max denominator z = ", sprintf("%.2f", max_den_z),
            "), so no percentage may be attached to any estimate.")
        else paste0(n_reportable, " of ", n_fam, " families pass the relative-magnitude ",
                    "gate; only those may be stated as a percentage."),
        paste0("A negative gradient is NOT evidence that exposure effects concentrate ",
               "in weakly controlled VMRs (variance budget, AGENTS.md 7.10).")),
        strwrap, width = 120)), collapse = "\n")

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
             caption = CAP_E) +
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
                      paste0("stage A FDR-surviving VMR x exposure pairs (",
                             paste(sprintf("%d %s", as.integer(table(factor(ever$region, levels = regions))),
                                           regions), collapse = " / "), ")"),
                      data_dir)
}

## ================================= GREML simulation benchmark (Module 02b)
##
## Recovery of ABSOLUTE local h2 by GCTA-GREML on SIMULATED phenotypes: arm 2
## on real AA cis-window genotypes at each region's design n, arm 1 the v1
## AR(1) design (out of regime). No observed methylation locus is estimated, so
## nothing here is a PVE for any VMR (AGENTS.md 3, 7.2). The panel shows what
## the ordering-vs-level distinction looks like for an estimator that is not
## the elastic net Module 02 retired.
##
## A production build requires the accepted 02b runs. A draft build
## (--out-dir) may read the -UNACCEPTED collation instead, and says so on the
## panel; it never reaches a sealed run that way.
GREML_MODULE <- "02b_greml_simulation_benchmark"
gcfg <- load_config("greml_benchmark")
greml_cells <- rbind(data.table(cohort = gcfg$arm1$run_cohort, region = gcfg$arm1$run_region),
                     data.table(cohort = gcfg$arm2$cohort, region = unlist(gcfg$arm2$regions)))
greml_ok <- all(vapply(seq_len(nrow(greml_cells)), function(i) !is.null(tryCatch(
    require_accepted_upstream(GREML_MODULE, greml_cells$cohort[i], greml_cells$region[i]),
    error = function(e) NULL)), logical(1)))
gcomb <- file.path(V2_ROOT, GREML_MODULE, "_m", "combined")
gsuffix <- if (greml_ok) "" else "-UNACCEPTED"
if (!greml_ok && is.null(opts$out_dir)) {
    stop("02b_greml_simulation_benchmark has no accepted run for every cell; ",
         "the GREML supplement cannot enter a sealed figure run.")
}
gm_f <- file.path(gcomb, paste0("greml-benchmark-recovery-metrics", gsuffix, ".tsv"))
if (file.exists(gm_f)) {
    gm <- fread(gm_f)
    gsp <- fread(file.path(gcomb, paste0("greml-benchmark-spearman", gsuffix, ".tsv")))
    stopifnot(all(gm$simulated_phenotypes_only %in% c(TRUE, "TRUE")),
              all(gm$absolute_pve_interpretation_allowed_for_observed_loci %in% c(FALSE, "FALSE")))
    PRIMARY_REML <- unique(gm[reml_role == "primary", reml_mode])
    stopifnot(length(PRIMARY_REML) == 1L)
    a2 <- gm[cell_type == "region_design_n" & reml_mode == PRIMARY_REML]
    a1 <- gm[cell_type == "simulated_n"]
    a2[, region := as_region(cell)]
    draft_note <- if (greml_ok) "" else "DRAFT: built from UNACCEPTED 02b runs. "

    pG1 <- ggplot(a2, aes(h2_nominal, mean_estimate, colour = region)) +
        geom_abline(slope = 1, intercept = 0, colour = PAL_NULL, linewidth = 0.35) +
        geom_line(linewidth = 0.45) + geom_point(size = 1.1) +
        facet_wrap(~ architecture, nrow = 1) +
        scale_colour_manual(values = REGION_COLORS, name = NULL) +
        labs(x = "Simulated local h2", y = "Mean REML estimate") +
        BASE_THEME + NO_TITLES + GRID_Y +
        theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0))
    pG2 <- ggplot(a2, aes(h2_nominal, coverage95, colour = region)) +
        geom_hline(yintercept = 0.95, colour = PAL_NULL, linewidth = 0.35, linetype = 2) +
        geom_line(linewidth = 0.45) + geom_point(size = 1.1) +
        facet_wrap(~ architecture, nrow = 1) +
        scale_colour_manual(values = REGION_COLORS, guide = "none") +
        scale_y_continuous(labels = percent_format(accuracy = 1)) +
        labs(x = "Simulated local h2", y = "95% CI coverage") +
        BASE_THEME + NO_TITLES + GRID_Y
    pG3 <- if (nrow(a1)) ggplot(a1, aes(as.numeric(cell), bias, colour = reml_mode)) +
        geom_hline(yintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
        geom_line(linewidth = 0.45) + geom_point(size = 1.1) +
        scale_x_log10() +
        scale_colour_manual(values = c(PAL_BLUE, PAL_TAN), name = NULL) +
        labs(x = "Simulated sample size (log scale)", y = "Mean bias") +
        BASE_THEME + NO_TITLES + GRID_Y +
        theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0))
    else NULL
    CAP_G <- paste(strwrap(paste0(
        draft_note, "Simulated phenotypes only; no observed locus is estimated. ",
        "a-b: real AA cis-window genotypes at each region's design n, primary mode ",
        PRIMARY_REML, ".", if (!is.null(pG3)) paste0(
            " c: the v1 AR(1) design, out of regime for real cis-windows; it does ",
            "not transfer to the cohort.") else ""), width = 115), collapse = "\n")
    cap_theme <- theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"),
                       plot.caption.position = "plot")
    ## The caption rides the last panel, whichever that is.
    if (is.null(pG3)) pG2 <- pG2 + labs(caption = CAP_G) + cap_theme
    else pG3 <- pG3 + labs(caption = CAP_G) + cap_theme
    SG <- paste0("figureS_greml_benchmark", arm)
    gfig <- if (is.null(pG3)) (pG1 / pG2) else (pG1 / pG2 / pG3)
    save_figure(gfig + fig_tags() & TAG_THEME, SG, width = FIG_WIDTH_FULL,
                height = if (is.null(pG3)) 5.0 else 7.4, fig_dir = fig_dir)
    gruns <- unique(gm$source_run_id)
    write_source_data(a2, paste0(SG, "_panel_a"), gruns,
                      paste0("_m/combined/greml-benchmark-recovery-metrics", gsuffix, ".tsv"),
                      SCRIPT, paste0("arm 2 (observed_AA_cis), reml_mode == ", PRIMARY_REML,
                                     "; mean estimate vs simulated h2"), data_dir)
    write_source_data(a2, paste0(SG, "_panel_b"), gruns,
                      paste0("_m/combined/greml-benchmark-recovery-metrics", gsuffix, ".tsv"),
                      SCRIPT, paste0("arm 2, reml_mode == ", PRIMARY_REML,
                                     "; 95% CI coverage of the realized simulated h2"), data_dir)
    if (nrow(a1)) write_source_data(a1, paste0(SG, "_panel_c"), gruns,
                      paste0("_m/combined/greml-benchmark-recovery-metrics", gsuffix, ".tsv"),
                      SCRIPT, "arm 1 (ar1_out_of_regime), both REML modes; out of regime", data_dir)
    write_source_data(gsp, paste0(SG, "_spearman"), gruns,
                      paste0("_m/combined/greml-benchmark-spearman", gsuffix, ".tsv"),
                      SCRIPT, "Spearman(truth, estimate) with cluster-bootstrap CI, every cell and mode",
                      data_dir)
} else {
    message("[skip] figureS_greml_benchmark: no ", basename(gm_f))
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
