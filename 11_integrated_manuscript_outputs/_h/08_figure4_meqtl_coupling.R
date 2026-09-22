#### 11 / Figure 4: CpG meQTL burden and transcription/splicing coupling ####
##
## AGENTS.md 11 assigns Figure 4 "CpG meQTL burden and expression/splicing
## coupling". Two owning modules: 05_cpg_meqtl_burden (AGENTS.md 7.5) and
## 07_transcription_splicing_coupling (AGENTS.md 7.6).
##
## Panels, in rendered order
##   a  meQTL-positive CpG fraction across deciles of the local-control rank
##   b  the burden model coefficient per region, with the distal-null genomic
##      inflation lambda stated as QC on the same panel
##   c  the coupling tests: modality x predictor x region, with the number of
##      coupled VMRs printed beside each estimate
##   d  the gradient behind panel c -- coupled-VMR fraction across deciles
##
## WHAT THIS FIGURE MAY AND MAY NOT SAY
## Module 07's permitted claim is exactly: "genetically regulated VMRs are more
## frequently transcriptionally coupled". AGENTS.md 7.6 forbids the mediation
## reading -- nothing here says methylation mediates a genetic effect on
## expression or splicing, and no axis, label or caption may imply it.
## Module 05 is CONVERGENT evidence, not independent replication: the meQTLs
## are mapped in the same donors as the score (AGENTS.md 7.5).
##
## PSI is not uniform and the panel must not let it look uniform: strong in
## caudate (227 coupled VMRs), thin in DLPFC (24) and null in hippocampus (7).
## `n_coupled` is therefore printed, not hidden in a caption.
##
## The ABC modality is underpowered (312 tested VMRs, 19 coupled, q 0.07-0.32)
## and is NOT a claim. It is excluded from the panel and kept in source data.
##
## Usage:
##   Rscript 08_figure4_meqtl_coupling.R --cohort AA --run-id fig-all-YYYYMMDD
##   Rscript 08_figure4_meqtl_coupling.R --cohort AA --run-id smoke \
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

SCRIPT <- "11_integrated_manuscript_outputs/_h/08_figure4_meqtl_coupling.R"

## ------------------------------------------------------------- upstream
CMB <- vapply(regions, function(r)
    require_accepted_upstream("05_cpg_meqtl_burden", cohort, r)$run_id, character(1))
TSC <- vapply(regions, function(r)
    require_accepted_upstream("07_transcription_splicing_coupling", cohort, r)$run_id,
    character(1))
cmb_dir <- function(r) file.path(V2_ROOT, "05_cpg_meqtl_burden", "_m", "runs",
                                 CMB[[r]], "results")
tsc_dir <- function(r) file.path(V2_ROOT, "07_transcription_splicing_coupling",
                                 "_m", "runs", TSC[[r]], "results")

by_region <- function(fun) rbindlist(lapply(regions, function(r) {
    d <- fun(r); if (is.null(d) || nrow(d) == 0) return(NULL); d[, region := r][]
}), fill = TRUE)

burden  <- by_region(function(r) fread(file.path(cmb_dir(r), "vmr-meqtl-burden.tsv")))
bmodel  <- by_region(function(r) fread(file.path(cmb_dir(r), "burden-primary-model.tsv")))
bdec    <- by_region(function(r) fread(file.path(cmb_dir(r), "burden-decision.tsv")))
coupling <- by_region(function(r) fread(file.path(tsc_dir(r), "coupling-tests.tsv")))
cdec    <- by_region(function(r) fread(file.path(tsc_dir(r), "coupling-decision.tsv")))

## -------------------------------------------------------------- guards
for (nm in c("burden", "bmodel", "coupling")) {
    banned <- intersect(c("h2_en_calibrated", "r_squared_cv", "h2_unscaled"),
                        names(get(nm)))
    if (length(banned) > 0) {
        stop("Retired quantity in ", nm, ": ", paste(banned, collapse = ", "),
             " (AGENTS.md 3).")
    }
}
if (!all(bdec$decision == "PASS_CPG_MEQTL_BURDEN_QC")) {
    stop("A Module 05 run does not carry PASS_CPG_MEQTL_BURDEN_QC.")
}
## The claim is about coupling FREQUENCY, never mediation. If Module 07 ever
## emits a mediation-flavoured term this figure must stop rather than render it.
if (any(grepl("mediat", names(coupling), ignore.case = TRUE)) ||
    any(grepl("mediat", coupling$modality, ignore.case = TRUE))) {
    stop("Module 07 emitted a mediation term; AGENTS.md 7.6 forbids that claim.")
}

for (d in list(burden, bmodel, coupling)) d[, region := as_region(region)]

## ---------------------- a. meQTL-positive CpG fraction across the rank
##
## Deciles of the rank, not a threshold (AGENTS.md 7.2). The y axis is the
## fraction of a VMR's TESTED CpGs with a significant cis meQTL, so the
## denominator is the tested universe and is reported with it.
burden[, decile := cut(local_snp_contribution_score, breaks = seq(0, 1, 0.1),
                       labels = 1:10, include.lowest = TRUE)]
bd <- burden[!is.na(decile) & n_tested_cpgs > 0,
             .(frac = sum(n_cpgs_with_sig_meqtl) / sum(n_tested_cpgs),
               n_vmrs = .N, n_tested_cpgs = sum(n_tested_cpgs)),
             by = .(region, decile)]

pA <- ggplot(bd, aes(as.integer(decile), frac, colour = region)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.3) +
    scale_colour_manual(values = REGION_COLORS, name = NULL) +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = "Decile of local SNP contribution rank",
         y = "Tested CpGs with\na significant cis meQTL") +
    BASE_THEME + NO_TITLES +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0))

## -------------------------------- b. the burden model, and its inflation QC
bm <- bmodel[term == "local_snp_contribution_score_z"]
stopifnot(nrow(bm) == length(regions))
bm[, `:=`(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)]
## Denominator and inflation ride INSIDE the panel, anchored to the left edge.
## Anchored to the estimate they overflowed the half-width column, and a prose
## caption did the same -- at 6.5 pt in ~3.2 in only about 55 characters fit,
## and ggplot does not wrap a caption. The "convergent evidence, not
## independent replication" qualifier is in the panel's source data and belongs
## in the figure legend, which has room for a sentence.
bm[, note := paste0("n = ", label_comma()(n_vmrs),
                    "   \u03bb = ", sprintf("%.3f", genomic_inflation_lambda))]

pB <- ggplot(bm, aes(estimate, region, colour = region)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    errorbar_h(aes(xmin = lo, xmax = hi), linewidth = 0.45) +
    geom_point(size = 1.8) +
    geom_text(aes(label = note), x = -Inf, hjust = -0.06, nudge_y = 0.3,
              size = 2.2, colour = "grey35", show.legend = FALSE) +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_y_discrete(limits = rev, expand = expansion(add = c(0.5, 0.75))) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.10))) +
    labs(x = "meQTL burden per SD of rank\n(quasibinomial)", y = NULL,
         caption = "\u03bb = distal-null genomic inflation") +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(plot.caption = element_text(size = 6.5, hjust = 0, colour = "grey35"))

## ------------------------------------------- c. transcription/splicing coupling
MODALITY_LABELS <- c(expression_nearest_gene = "Expression\n(nearest gene)",
                     psi                     = "Splicing (PSI)")
PREDICTOR_LABELS <- c(local_genetic_control = "Local genetic control",
                      meqtl_proportion      = "meQTL proportion",
                      any_meqtl_support     = "Any meQTL support")

## ABC is underpowered and is not a claim (AGENTS.md 7.6): kept in source data,
## kept out of the panel.
cp <- coupling[modality %in% names(MODALITY_LABELS)]
cp[, `:=`(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se,
          mod = factor(MODALITY_LABELS[modality], levels = unname(MODALITY_LABELS)),
          pred = factor(PREDICTOR_LABELS[predictor],
                        levels = rev(unname(PREDICTOR_LABELS))))]
cp[, stars := sig_stars(q)]

pC <- ggplot(cp, aes(estimate, pred, colour = region)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    errorbar_h(aes(xmin = lo, xmax = hi),
               position = position_dodge(width = 0.68), linewidth = 0.42) +
    geom_point(size = 1.5, position = position_dodge(width = 0.68)) +
    geom_text(aes(label = stars), position = position_dodge(width = 0.68),
              hjust = -0.35, vjust = 0.75, size = 2.6, show.legend = FALSE) +
    facet_wrap(~ mod, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    labs(x = "Coupling log-odds per SD of predictor", y = NULL) +
    BASE_THEME + NO_TITLES + GRID_Y

## ------------------------------- d. how many VMRs each estimate rests on
##
## The denominators, on the same figure as the estimates. PSI is strong in
## caudate, thin in DLPFC and null in hippocampus, and a forest plot alone
## hides that.
den <- unique(cp[, .(region, mod, n, n_coupled)])
den[, frac := n_coupled / n]

pD <- ggplot(den, aes(region, n_coupled, fill = region)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = paste0(label_comma()(n_coupled), "\n",
                                 percent(frac, accuracy = 0.1), " of ",
                                 label_comma()(n))),
              vjust = -0.25, size = 2.2, colour = "black", lineheight = 0.95) +
    facet_wrap(~ mod, nrow = 1) +
    scale_fill_manual(values = REGION_COLORS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
    labs(x = NULL, y = "Coupled VMRs") +
    BASE_THEME + NO_TITLES +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

## ------------------------------------------------------------------ assemble
arm <- if (cohort == "AA") "" else paste0("_", cohort)
STEM <- paste0("figure4_meqtl_burden_coupling", arm)

figure <- ((pA | pB) + plot_layout(widths = c(1, 0.8))) / pC / pD +
    plot_layout(heights = c(1.0, 0.95, 0.9)) +
    fig_tags() & TAG_THEME

save_figure(figure, STEM, width = FIG_WIDTH_FULL, height = 7.6,
            fig_dir = fig_dir)

## ---------------------------------------------------------- source data
runs_used <- c(unname(CMB), unname(TSC))
sd <- function(dt, nm, tbl, filt) {
    write_source_data(dt, paste0(STEM, "_", nm), runs_used, tbl, SCRIPT,
                      filt, data_dir)
}
sd(bd, "panel_a", "05: results/vmr-meqtl-burden.tsv",
   "n_tested_cpgs > 0; deciles of local_snp_contribution_score; fraction is summed sig CpGs over summed TESTED CpGs")
sd(bm[, .(region, model, term, estimate, se, z, p, n_vmrs, dispersion,
          genomic_inflation_lambda)],
   "panel_b", "05: results/burden-primary-model.tsv",
   "term == 'local_snp_contribution_score_z'; convergent evidence, not independent replication")
sd(coupling[, .(region, modality, predictor, n, n_coupled, estimate, se, z, p, q,
                covariates)],
   "panel_c", "07: results/coupling-tests.tsv",
   "all modalities INCLUDING expression_abc, which is underpowered (n=312, 19 coupled) and is excluded from the rendered panel and is not a claim")
sd(den, "panel_d", "07: results/coupling-tests.tsv",
   "coupled-VMR counts and tested universe per modality x region")
sd(cdec, "coupling_decision", "07: results/coupling-decision.tsv",
   "the accepted Module 07 decision row")
sd(bdec, "burden_decision", "05: results/burden-decision.tsv",
   "the accepted Module 05 decision row")

message("[done] Figure 4 written to ", fig_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
