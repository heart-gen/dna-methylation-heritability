#### 11_integrated_manuscript_outputs / 00_nongenetic_panels: shared panels ####
##
## Three panels are rendered in two figures each. The age gradient and the
## exposure gradients open Figure 6, where the low end of the local
## genetic-control axis is read against non-genetic variation and then against
## complex-trait loci, and they also anchor their own supplementary figures
## (figureS_aging_axis a, figureS_environmental_axis), which keep the
## sensitivity detail. The psychiatric-vs-other contrast left Figure 6 for
## figureS_schizophrenia_application d on 2026-10-22.
##
## One builder per panel, sourced by 09_figure6_gwas_nongenetic_axis.R and
## 11_supplementary_figures.R. Two copies of a panel is how v1's figures
## drifted apart (00_figure_theme.R, header), so the loaders and the guards
## that used to sit inline in the supplementary builder live here once.
##
## Not a standalone script. load_*() resolves an upstream through the README
## of record and returns its tables with the guards applied; panel_*() takes
## those tables and returns the ggplot and the legend sentence derived from
## them. Nothing here writes a file: figure assembly, legend_note() and
## write_source_data() stay with the builder that owns each figure.

.by_region <- function(regions, fun) rbindlist(lapply(regions, function(r) {
    d <- fun(r); if (is.null(d) || nrow(d) == 0) return(NULL); d[, region := r][]
}), fill = TRUE)

## Stage 17's primary model, named once for the two builders that filter on it.
STAGE17_PRIMARY_MODEL <- "A_technical"

## ===================================================== Aging (Module 09b)
##
## The VERDICT -- the cross-region token, which regions are supported, which
## arms fail -- is read out of the run's own decision tables and never written
## into a string here. A caption that states a per-region verdict as prose goes
## stale silently the next time the verdict changes, and this one did while the
## module 09b scMD-gate correction was in flight. AGENTS.md 7.11 wants a
## panel's numbers traceable to a run, table, script and filter.
load_aging <- function(cohort, regions) {
    AGE <- vapply(regions, function(r)
        require_accepted_upstream("09b_aging_application", cohort, r)$run_id,
        character(1))
    age_dir <- function(r) file.path(V2_ROOT, "09b_aging_application", "_m",
                                     "runs", AGE[[r]], "results")
    per <- fread(file.path(V2_ROOT, "09b_aging_application", "_m", "combined",
                           sprintf("aging-per-region-%s.tsv", cohort)))
    gat <- .by_region(regions, function(r)
        fread(file.path(age_dir(r), "gating-sensitivities.tsv")))
    qrt <- .by_region(regions, function(r)
        fread(file.path(age_dir(r), "axis-quartile-summary.tsv")))

    ## The stage-05 cross-region decision row, which is where the token lives.
    xr <- fread(file.path(V2_ROOT, "09b_aging_application", "_m", "combined",
                          sprintf("aging-cross-region-decision-%s.tsv", cohort)))
    if (nrow(xr) != 1L) {
        stop("Expected one aging cross-region decision row for ", cohort,
             ", found ", nrow(xr))
    }
    token <- as.character(xr$aging_axis_association)
    if (!nzchar(token) || is.na(token)) {
        stop("Aging cross-region decision carries no aging_axis_association token.")
    }

    stopifnot(all(per$cross_sectional_design %in% c(TRUE, "TRUE")))
    stopifnot(all(per$causal_interpretation_allowed %in% c(FALSE, "FALSE")))
    ## The 2026-10-03 magnitude gate (T28 decision B): a region whose family
    ## mean is not separated from zero reports a SIGNED TEST only, never a
    ## proportional change. The label says which, from the run's own column, so
    ## a reader cannot take the point estimate of a gated region as a percentage.
    per[, mag_ok := primary_relative_magnitude_reportable %in% c(TRUE, "TRUE")]
    per[, p_lab := paste0("P = ", signif(primary_p, 2),
                          fifelse(mag_ok, "", " · signed test only"))]
    per[, region := as_region(region)]
    gat[, region := as_region(region)]
    list(runs = unname(AGE), per = per, gat = gat, qrt = qrt, token = token)
}

#' The primary age gradient per region: Figure 6a and figureS_aging_axis a.
panel_aging_primary <- function(per) {
    plot <- ggplot(per, aes(primary_estimate, region, colour = region)) +
        geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
        errorbar_h(aes(xmin = primary_ci_lower, xmax = primary_ci_upper),
                   linewidth = 0.45) +
        geom_point(size = 1.8) +
        geom_text(aes(label = p_lab),
                  x = Inf, hjust = 1.05, vjust = -1.0, size = 2.2,
                  colour = "grey35") +
        scale_colour_manual(values = REGION_COLORS, guide = "none") +
        scale_y_discrete(limits = rev, expand = expansion(add = c(0.5, 0.8))) +
        labs(x = "Age-effect gradient per SD of rank", y = NULL) +
        BASE_THEME + NO_TITLES + GRID_Y
    note <- paste0(
        "Gradient in the debiased squared age effect along the local SNP ",
        "contribution rank, ratio scale, with 95% intervals. Cross-sectional design: ",
        "age-associated differences, never change with age. Signed test only: the ",
        "magnitude gate withholds a proportional reading (", sum(!per$mag_ok), " of ",
        nrow(per), " regions).")
    list(plot = plot, note = note)
}

## =============================================== Environmental (Module 10)
##
## Exploratory and gated on acceptance (AGENTS.md 6). The supplementary builder
## asks to SKIP when the acceptance rows are missing, so one blocked upstream
## does not stop the other supplements; Figure 6 cannot be assembled without
## the panel and stops instead.
load_environmental <- function(cohort, regions, run_dir,
                               skip_if_unaccepted = FALSE) {
    runs <- tryCatch(
        vapply(regions, function(r)
            require_accepted_upstream("10_environmental_exploratory", cohort, r)$run_id,
            character(1)),
        error = function(e) {
            if (!skip_if_unaccepted) stop(e)
            message("[skip] figureS_environmental_axis: ", conditionMessage(e))
            NULL
        })
    if (is.null(runs)) return(list(ok = FALSE))

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
    ## Readable family labels: "nicotine @ all" was a code identifier.
    STRATUM_LABELS <- c(all = "all donors", schizophrenia = "schizophrenia cases")
    eax[, fam := paste0(toupper(substr(gsub("_", " ", exposure), 1, 1)),
                        substring(gsub("_", " ", exposure), 2), " (",
                        fifelse(stratum %in% names(STRATUM_LABELS),
                                STRATUM_LABELS[stratum], gsub("_", " ", stratum)),
                        ")")]
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
    cap <- paste(
        paste0("Exploratory; x axis clipped to [", LIM[1], ", ", LIM[2],
               "], estimates outside it printed in full. Stars give the FDR call; ",
               "? marks a call that is undetermined across bootstrap seeds."),
        if (n_reportable == 0L) paste0(
            "Signed gradients only: no family passes the relative-magnitude gate (0 of ",
            n_fam, "; max denominator z = ", sprintf("%.2f", max_den_z),
            "), so no percentage may be attached to any estimate.")
        else paste0(n_reportable, " of ", n_fam, " families pass the relative-magnitude ",
                    "gate; only those may be stated as a percentage."),
        paste0("Whether active, low-local-control VMRs respond more to non-genetic ",
               "influences cannot be distinguished from the larger residual variance ",
               "that remains when local genetic contribution is lower."))

    list(ok = TRUE, runs = unique(eax$run_id), eax = eax, ever = ever,
         rel = rel, LIM = LIM, cap = cap)
}

#' Stage B exposure gradients: Figure 6b and figureS_environmental_axis.
#'
#' `legend = FALSE` drops the region key: inside Figure 6 the key is drawn once
#' for the whole figure from panel e, and a second copy over this panel is what
#' patchwork produced when both were left on. Standalone, the panel keeps it.
panel_environmental <- function(rel, LIM, legend = TRUE) {
    ggplot(rel, aes(primary_beta, fam, colour = region)) +
        geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
        errorbar_h(aes(xmin = primary_ci_lower, xmax = primary_ci_upper),
                   position = position_dodge(width = 0.72), linewidth = 0.4,
                   show.legend = FALSE) +
        geom_point(size = 1.5, position = position_dodge(width = 0.72)) +
        geom_text(aes(x = pmin(primary_ci_upper, LIM[2]), label = stars),
                  position = position_dodge(width = 0.72),
                  hjust = -0.25, vjust = 0.75, size = 2.5, show.legend = FALSE) +
        geom_text(data = rel[offscale == TRUE], aes(label = off_lab),
                  x = LIM[1], hjust = -0.03, vjust = -0.9, size = 2.0,
                  show.legend = FALSE) +
        coord_cartesian(xlim = LIM) +
        scale_colour_manual(values = REGION_COLORS, name = NULL,
                            guide = if (legend) "legend" else "none") +
        ## Two lines: at Figure 6's half width the one-line title was clipped.
        labs(x = "Exposure-explained variance gradient\nper SD of rank",
             y = NULL, caption = paste0(SIG_KEY, "   ? undetermined")) +
        BASE_THEME + NO_TITLES + GRID_Y + KEY_THEME +
        theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0))
}

## ====================================== Psychiatric vs other traits (Module 09)
##
## Stage 17's category contrast: psychiatric-neurologic traits are
## indistinguishable from the rest. A null with one number per region, it was
## the main figure's panel c until the age and exposure gradients took its
## place; the headline it supports (schizophrenia is a typical trait) still
## stands in Figure 6d, and the contrast itself is figureS_schizophrenia_application d.
panel_psychiatric_vs_other <- function(traits, summ, model_id) {
    psy <- traits[in_distribution == TRUE & model == model_id & is.finite(estimate)]
    ## Short labels: "Psychiatric-neurologic" and "All other traits" overlapped at
    ## this panel's width and squeezed the facet strips into "ppocamp". The full
    ## grouping is in the caption and the source data.
    psy[, grp := fifelse(category == "Psychiatric-neurologic",
                         "Psychiatric", "Other")]
    psy[, grp := factor(grp, levels = c("Psychiatric", "Other"))]
    psy[, region := as_region(region)]
    wil <- summ[model == model_id, .(region, psychiatric_vs_other_wilcoxon_p)]
    wil[, region := as_region(region)]
    wil[, lab := paste0("P = ", signif(psychiatric_vs_other_wilcoxon_p, 2))]

    ## Region on x, group as fill, and NO facets. Faceting by region gave each
    ## panel about 1.2 in for two category labels and a strip, so the labels
    ## collided and "Hippocampus" was clipped to "ppocampu". One panel also puts
    ## the two groups side by side within a region, which is the comparison.
    plot <- ggplot(psy, aes(region, estimate, fill = grp)) +
        geom_hline(yintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
        geom_boxplot(width = 0.62, outlier.size = 0.4, linewidth = 0.3,
                     position = position_dodge(width = 0.72)) +
        geom_text(data = wil, aes(x = region, y = Inf, label = lab),
                  inherit.aes = FALSE, vjust = 1.5, size = 2.3, colour = "grey35") +
        scale_fill_manual(values = c(Psychiatric = PAL_ALT, Other = PAL_NULL),
                          name = NULL) +
        ## Headroom so the P label clears the topmost outlier.
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
        labs(x = NULL, y = "Axis estimate (score SD)") +
        BASE_THEME + NO_TITLES + GRID_Y +
        theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0))
    note <- paste(
        "Psychiatric, the psychiatric-neurologic GWAS category; Other, all remaining",
        "traits. P, two-sided Wilcoxon test of psychiatric against other traits.")
    list(plot = plot, psy = psy, wil = wil, note = note)
}
