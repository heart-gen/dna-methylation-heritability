#### 10 / Figure 2: relative local genetic control and held-out prediction ####
##
## AGENTS.md 11 assigns Figure 2 local SNP-explained variance and the secondary
## held-out prediction. Module 02's terminal decision is
## PASS_RELATIVE_GENETIC_CONTROL_FAIL_ABSOLUTE_LOCUS_PVE, so the *only*
## admissible endpoint is the relative score: absolute PVE percentages, any
## threshold (including 0.10), heritable/non-heritable groups, and the retired
## h2_en_calibrated are all prohibited. No main panel carries a PVE axis.
##
## What this figure can and cannot show
##   local_snp_contribution_score is a within-cell midrank percentile, so its
##   distribution is uniform BY CONSTRUCTION (observed sd 0.2887 = uniform).
##   Plotting that distribution, or its z-score, would show only the definition.
##   The evidence is instead in what the ranking agrees with:
##
##   a  the rank is reproduced by independent estimators, and is NOT explained
##      by locus geometry (the reviewer's first objection)
##   b  held-out local SNP prediction accuracy rises monotonically across the
##      rank -- the secondary endpoint of AGENTS.md 11
##
## Which prediction number panel b may carry
##   AGENTS.md 7.3 names ONE primary v2 prediction endpoint: `r2_pred_oof`,
##   the END-TO-END out-of-fold R-squared emitted by module 03, in which the
##   locus screen and the residualization are also learned inside the outer
##   training donors. Module 02 emits its own `r2_oof` from the nested CV
##   inside its joint-feature elastic net; that is a MODEL-LEVEL out-of-fold
##   statistic, and AGENTS.md 4 separates the two standards explicitly.
##
##   Until 2026-09-23 this panel plotted module 02's `r2_oof` under the axis
##   label "Held-out R2", which presented the weaker standard where the
##   manuscript claims the stronger, and left module 03 cited by no panel of
##   any figure. Panel b now reads `r2_pred_oof` from the accepted module 03
##   runs; module 02's `r2_oof` stays in panel a, relabelled as the
##   model-level statistic it is.
##   c  the rank is concordant across brain regions in VMRs called in both
##   d  genic context across the rank -- descriptive proportions only, no
##      enrichment test, no threshold. Repeat and repressive-chromatin
##      enrichment is Figure 3 / Module 04, whose enrichment model is a locked
##      PI decision (AGENTS.md 12); nothing here anticipates it.
##
## Denominators and exclusions ship as their own supplementary figure. AGENTS.md
## 11 requires they be reported, not that they occupy a main panel, and the
## exclusions/denominator table in 10_manuscript_tables.R carries the same
## numbers. Moving them out is what brings the main figure under a page.
##
## Usage:
##   Rscript 02_figure2_local_control.R --cohort AA --run-id fig-all-20260826

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(V2_ROOT, "11_integrated_manuscript_outputs", "_h",
                 "00_figure_theme.R"))

suppressPackageStartupMessages({
    library(patchwork)
    library(scales)
    library(GenomicRanges)
})

opts <- parse_v2_args(require = c("cohort", "run_id"))
cohort <- opts$cohort
regions <- load_config("cohorts")$regions

module_root <- file.path(V2_ROOT, "11_integrated_manuscript_outputs")
## --out-dir renders a review draft outside the immutable run tree, so a panel
## can be iterated on without minting and sealing a run. Same flag as
## 06_figure_region_donor_generalization.R.
run_dir  <- if (!is.null(opts$out_dir)) opts$out_dir else
    file.path(module_root, "_m", "runs", opts$run_id)
fig_dir  <- file.path(run_dir, "figures")
data_dir <- file.path(run_dir, "source_data")

## Resolved through the acceptance gate (AGENTS.md 6), never from a run-ID
## template. The template this replaced read `lgv-{cohort}-{region}-20260823`,
## which Module 02 retired on 2026-09-17 in favour of the rescored
## `lgv-AA-*-rescore-20260913` runs -- and because nothing checked, the last
## build of this figure (fig-all-20260826-a) shipped on the retired score and
## could not be cited.
LGV_RUN <- local({
    cache <- new.env(parent = emptyenv())
    function(r) {
        if (is.null(cache[[r]])) {
            cache[[r]] <- require_accepted_upstream("02_local_genetic_variance",
                                                    cohort = cohort,
                                                    region = r)$run_id
        }
        cache[[r]]
    }
})
lgv_file <- function(r) file.path(
    V2_ROOT, "02_local_genetic_variance", "_m", "runs", LGV_RUN(r),
    "results", "combined",
    sprintf("local-genetic-control-%s-%s-vmrs.tsv", cohort, r))

all_rows <- rbindlist(lapply(regions, function(r) {
    d <- fread(lgv_file(r)); d[, region := r][]
}), fill = TRUE)

## Guard the prohibition at runtime rather than trusting the script to stay
## clean: a retired column reappearing upstream must stop the figure, not
## silently enter it.
banned <- intersect(c("h2_en_calibrated", "positive_signal"), names(all_rows))
if (length(banned) > 0) {
    stop("Retired quantity present in the module 02 contract: ",
         paste(banned, collapse = ", "), " (AGENTS.md 3).")
}
stopifnot(all(all_rows$absolute_pve_interpretation_allowed == FALSE))
stopifnot(all(all_rows$local_snp_contribution_score_basis == "pve_cis_joint_unbounded"))

elig <- all_rows[local_genetic_control_eligible == TRUE]
elig[, region := as_region(region)]

## -------------------- module 03: the end-to-end out-of-fold prediction endpoint
##
## Resolved through the same acceptance gate as module 02 (AGENTS.md 6), so the
## run ID lands in this figure's source data and from there in the run manifest.
## The join is on vmr_id and is safe on identity grounds only because both
## modules key on the same accepted catalog: assert the vmr_set_id agrees rather
## than trusting it. Module 03's accepted runs consumed the PRE-rescore module 02
## run for their locus screen, which is why the check is on vmr_set_id -- the
## catalog -- and not on the upstream module 02 run ID. `r2_pred_oof` is a
## genotype-to-phenotype quantity and carries no score in it, so the rescore
## does not touch it.
LSP_RUN <- local({
    cache <- new.env(parent = emptyenv())
    function(r) {
        if (is.null(cache[[r]])) {
            cache[[r]] <- require_accepted_upstream("03_local_snp_prediction",
                                                    cohort = cohort,
                                                    region = r)$run_id
        }
        cache[[r]]
    }
})
lsp_file <- function(r) file.path(
    V2_ROOT, "03_local_snp_prediction", "_m", "runs", LSP_RUN(r),
    "results", "combined",
    sprintf("oof-prediction-%s-%s-vmrs.tsv", cohort, r))

pred <- rbindlist(lapply(regions, function(r) {
    d <- fread(lsp_file(r))
    legacy <- intersect(c("r_squared_cv", "h2_unscaled", "h2_en_calibrated"),
                        names(d))
    if (length(legacy) > 0) {
        stop("Module 03 table for ", r, " carries retired metric(s): ",
             paste(legacy, collapse = ", "), " (AGENTS.md 3).")
    }
    if (!"r2_pred_oof" %in% names(d)) {
        stop("Module 03 table for ", r, " has no r2_pred_oof column; AGENTS.md ",
             "7.3 names it as the primary prediction endpoint.")
    }
    want <- unique(all_rows[region == r]$vmr_set_id)
    got <- unique(as.character(d$vmr_set_id))
    if (!identical(sort(want), sort(got))) {
        stop("vmr_set_id disagrees between module 02 (", paste(want, collapse = ","),
             ") and module 03 (", paste(got, collapse = ","), ") for ", r,
             "; the panel b join would cross VMR catalogs.")
    }
    d[, .(vmr_id, region = r, r2_pred_oof)]
}))
pred[, region := as_region(region)]

n_before <- nrow(elig)
elig <- merge(elig, pred, by = c("vmr_id", "region"), all.x = TRUE)
if (nrow(elig) != n_before) {
    stop("Module 03 join changed the eligible row count (", n_before, " -> ",
         nrow(elig), "); vmr_id is not unique within region.")
}
message("[join] end-to-end OOF prediction matched ",
        sum(!is.na(elig$r2_pred_oof)), " of ", nrow(elig), " eligible VMRs")

SCRIPT <- "11_integrated_manuscript_outputs/_h/02_figure2_local_control.R"
runs_used <- vapply(regions, LGV_RUN, "")
runs_used_pred <- c(runs_used, vapply(regions, LSP_RUN, ""))
FILTER <- "local_genetic_control_eligible == TRUE"

## ------------------------- A. what the rank agrees with, and what it does not
##
## Two contrasting groups on one axis. Independent estimators of local genetic
## control should track the rank; locus geometry should not, because a ranking
## driven by SNP count or LD would be an artifact rather than a signal.
## The two `*_oof` entries are module 02's own nested-CV statistics, so they are
## named for the standard they meet (AGENTS.md 4: model-level, not end-to-end).
## The bare label "Held-out R2" conflated them with panel b's endpoint.
CONCORD <- c(bslmm_pve = "BSLMM PVE",
             he_h2     = "Haseman-Elston",
             rho2_oof  = "Model-level OOF \u03c1\u00b2",
             r2_oof    = "Model-level OOF R\u00b2")
GEOMETRY <- c(num_snps  = "cis SNPs",
              p_eff     = "Effective rank",
              ld_metric = "LD")

spearman_ci <- function(x, y) {
    ok <- stats::complete.cases(x, y); x <- x[ok]; y <- y[ok]
    n <- length(x)
    rho <- stats::cor(x, y, method = "spearman")
    ## Fisher z on the Spearman rho, with the Bonett-Wright standard error.
    se <- sqrt((1 + rho^2 / 2) / (n - 3))
    z <- atanh(rho)
    list(rho = rho, lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se), n = n)
}

conc <- rbindlist(lapply(names(c(CONCORD, GEOMETRY)), function(v) {
    elig[, {
        s <- spearman_ci(local_snp_contribution_score, get(v))
        .(variable = v, rho = s$rho, lo = s$lo, hi = s$hi, n = s$n)
    }, by = region]
}))
conc[, group := fifelse(variable %in% names(CONCORD),
                        "Estimators", "Geometry")]
conc[, label := c(CONCORD, GEOMETRY)[variable]]
conc[, label := factor(label, levels = rev(c(CONCORD, GEOMETRY)))]
conc[, group := factor(group, levels = c("Estimators", "Geometry"))]

pA <- ggplot(conc, aes(rho, label, colour = region)) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    errorbar_h(aes(xmin = lo, xmax = hi),
               position = position_dodge(width = 0.62), linewidth = 0.42) +
    geom_point(size = 1.5, position = position_dodge(width = 0.62)) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = REGION_COLORS, name = NULL) +
    scale_x_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Spearman correlation with local SNP contribution rank", y = NULL) +
    BASE_THEME + NO_TITLES + GRID_Y +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          strip.text.y.right = element_text(angle = -90, face = "bold"),
          plot.margin = margin(5, 12, 5, 8))

## ------------------ B. held-out local SNP prediction across the rank
##
## The secondary endpoint, and it is module 03's `r2_pred_oof` (AGENTS.md 7.3).
## Deciles of the rank, not thresholds: no cut point is claimed, and the axis is
## prediction accuracy, never variance explained.
##
## AGENTS.md 7.3 also requires that negative `r2_pred_oof` be RETAINED and never
## swapped for `cor2_oof` when it is unfavourable. There is no floor and no
## drop here: median and quartiles are taken on the raw column, and the count of
## negative loci per decile ships in the panel's source data so the retention is
## auditable from the table rather than asserted in a comment. Most low-decile
## loci are negative, which is the honest reading of "not imputable".
dec <- copy(elig)
dec[, decile := cut(local_snp_contribution_score, breaks = seq(0, 1, 0.1),
                    labels = 1:10, include.lowest = TRUE)]
dec_sum <- dec[!is.na(decile), .(
    n = .N,
    n_r2_negative = sum(r2_pred_oof < 0, na.rm = TRUE),
    n_r2_missing = sum(is.na(r2_pred_oof)),
    median = median(r2_pred_oof, na.rm = TRUE),
    q25 = quantile(r2_pred_oof, 0.25, na.rm = TRUE),
    q75 = quantile(r2_pred_oof, 0.75, na.rm = TRUE)), by = .(region, decile)]

pB <- ggplot(dec_sum, aes(as.integer(decile), median, colour = region, fill = region)) +
    geom_hline(yintercept = 0, colour = PAL_NULL, linewidth = 0.35) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.3) +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_fill_manual(values = REGION_COLORS, guide = "none") +
    scale_x_continuous(breaks = 1:10) +
    labs(x = "Decile of local SNP contribution rank",
         y = expression(atop("Held-out local SNP prediction",
                             R^2 ~ "(end-to-end out-of-fold)"))) +
    BASE_THEME + NO_TITLES

## ------------------------------------ C. cross-region rank concordance
##
## VMR sets are called per region, so loci are matched by genomic overlap
## rather than by vmr_set_id, which is region-specific by design.
pair_concordance <- function(r1, r2) {
    a <- elig[region == REGION_LABELS[[r1]]]
    b <- elig[region == REGION_LABELS[[r2]]]
    ga <- GRanges(a$chrom, IRanges(a$start, a$end))
    gb <- GRanges(b$chrom, IRanges(b$start, b$end))
    ov <- findOverlaps(ga, gb)
    ## One VMR can overlap several in the other region; keep the widest overlap
    ## so each locus contributes once.
    w <- width(pintersect(ga[queryHits(ov)], gb[subjectHits(ov)]))
    dt <- data.table(qi = queryHits(ov), si = subjectHits(ov), w = w)
    setorder(dt, qi, -w)
    dt <- unique(dt, by = "qi")
    setorder(dt, si, -w)
    dt <- unique(dt, by = "si")
    s <- spearman_ci(a$local_snp_contribution_score[dt$qi],
                     b$local_snp_contribution_score[dt$si])
    data.table(pair = paste(REGION_LABELS[[r1]], "vs", REGION_LABELS[[r2]]),
               rho = s$rho, lo = s$lo, hi = s$hi, n = s$n)
}
pairs_dt <- rbindlist(list(
    pair_concordance("caudate", "dlpfc"),
    pair_concordance("caudate", "hippocampus"),
    pair_concordance("dlpfc", "hippocampus")))
pairs_dt[, pair := factor(pair, levels = rev(pair))]

pC <- ggplot(pairs_dt, aes(rho, pair)) +
    errorbar_h(aes(xmin = lo, xmax = hi),
               colour = PAL_CHARCOAL, linewidth = 0.45) +
    geom_point(size = 1.8, colour = PAL_RUST) +
    geom_text(aes(x = 0.02, label = paste0("n = ", label_comma()(n))),
              hjust = 0, nudge_y = 0.28, size = 2.4, colour = "grey35") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Rank concordance in shared VMRs", y = NULL) +
    BASE_THEME + NO_TITLES

## ------------------------------- D. genic context across the rank
##
## Joined to the Module 01 QC run on vmr_id: both modules key on the same
## accepted catalog, so the join is exact rather than positional.
## Module 01's QC refresh: no acceptance row of its own, so it is named in
## 00_figure_theme.R rather than duplicated here and in Figure 1.
CTX_RUN <- function(r) QC_REFRESH_RUN(cohort, r)
ctx <- rbindlist(lapply(regions, function(r) {
    d <- fread(file.path(V2_ROOT, "01_vmr_catalog", "_m", "runs", CTX_RUN(r),
                         "qc", "distance_to_nearest_gene.tsv"))
    d[, region := r][, .(vmr_id, region, genomic_context)]
}))
ctx[, region := as_region(region)]

dec_ctx <- merge(dec[!is.na(decile), .(vmr_id, region, decile)], ctx,
                 by = c("vmr_id", "region"))
if (nrow(dec_ctx) == 0) stop("Genic context join produced no rows; check vmr_id keys.")
message("[join] genic context matched ", nrow(dec_ctx), " of ", nrow(dec), " scored VMRs")

CTX_LEVELS <- c("Promoter", "5\' UTR", "Exon", "Intron", "3\' UTR", "Intergenic")
ctx_sum <- dec_ctx[, .N, by = .(region, decile, genomic_context)]
ctx_sum[, frac := N / sum(N), by = .(region, decile)]
## Show the compartments carrying the signal; the rare UTR classes are in the
## source-data table rather than as near-zero lines nobody can read.
SHOW <- c("Promoter", "Intron", "Intergenic")
ctx_show <- ctx_sum[genomic_context %in% SHOW]
ctx_show[, genomic_context := factor(genomic_context, levels = SHOW)]

pCtx <- ggplot(ctx_show, aes(as.integer(decile), frac, colour = region)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.1) +
    facet_wrap(~ genomic_context, nrow = 1) +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_x_continuous(breaks = c(1, 5, 10)) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = "Decile of local SNP contribution rank", y = "VMRs") +
    BASE_THEME + NO_TITLES

## ------------------------------------------ E. denominators and exclusions
excl <- all_rows[, .(
    Eligible = sum(local_genetic_control_eligible),
    Excluded = sum(!local_genetic_control_eligible)), by = region]
excl[, region := as_region(region)]
excl_long <- melt(excl, id.vars = "region", variable.name = "status",
                  value.name = "n")

reasons <- all_rows[local_genetic_control_eligible == FALSE,
                    .N, by = local_genetic_control_exclusion_reason]
setorder(reasons, -N)

pDenom <- ggplot(excl_long, aes(region, n, fill = status)) +
    geom_col(width = 0.68) +
    geom_text(data = excl_long[, .(n = sum(n),
                                   lab = paste0(label_comma()(n[status == "Eligible"]),
                                                " (", n[status == "Excluded"],
                                                " excluded)")), by = region],
              aes(label = lab, fill = NULL), vjust = -0.4, size = 2.4,
              colour = "black") +
    scale_fill_manual(values = c(Eligible = PAL_CHARCOAL, Excluded = PAL_NULL),
                      name = NULL) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()),
                       expand = expansion(mult = c(0, 0.14))) +
    labs(x = NULL, y = "VMRs tested") +
    BASE_THEME + NO_TITLES +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          axis.text.x = element_text(angle = 35, hjust = 1))

## ------------------------------------------------------------------ assemble
arm <- if (cohort == "AA") "" else paste0("_", cohort)
STEM <- paste0("figure2_local_genetic_control", arm)
figure <- (pA / (pB | pC) / pCtx) +
    plot_layout(heights = c(1.25, 1.0, 0.85)) +
    fig_tags() & TAG_THEME

save_figure(figure, STEM, width = FIG_WIDTH_FULL, height = 8.2,
            fig_dir = fig_dir)

## ------------------------------- supplement: denominators and exclusions
DENOM_STEM <- paste0("figureS_local_control_denominators", arm)
save_figure(pDenom + fig_tags() & TAG_THEME, DENOM_STEM,
            width = FIG_WIDTH_THREEQ, height = 3.4, fig_dir = fig_dir)

## --------------------------------------------------- supplement: audit only
##
## The unbounded joint estimate is the score's basis. It is shown ONLY as a
## diagnostic: module 02 failed its absolute-PVE gate, so no value here is
## interpretable as variance explained, and nothing downstream may use it.
pS <- ggplot(elig, aes(pve_cis_joint_unbounded, colour = region)) +
    geom_density(linewidth = 0.5) +
    geom_vline(xintercept = 0, colour = PAL_NULL, linewidth = 0.35,
               linetype = "dashed") +
    scale_colour_manual(values = REGION_COLORS, name = NULL) +
    labs(x = paste("Unbounded joint estimate (audit only; NOT interpretable as",
                   "\nabsolute variance explained -- module 02 failed its absolute-PVE gate)"),
         y = "Density") +
    BASE_THEME + NO_TITLES + theme(legend.position = "top")

save_figure(pS, paste0("figureS_local_control_audit_unbounded", arm),
            width = FIG_WIDTH_THREEQ, height = 3.4, fig_dir = fig_dir)

## ---------------------------------------------------------- source data
## Panel names follow the RENDERED tag: pA -> a, pB -> b, pC -> c, pCtx -> d.
sd <- function(dt, nm, tbl, filt, runs = runs_used) {
    write_source_data(dt, paste0(STEM, "_", nm),
                      runs, tbl, SCRIPT, filt, data_dir)
}
sd_denom <- function(dt, nm, tbl, filt) {
    write_source_data(dt, paste0(DENOM_STEM, "_", nm),
                      runs_used, tbl, SCRIPT, filt, data_dir)
}
## The audit panel ships as its own supplemental figure, so its source data is
## named for that figure rather than for Figure 2.
sd_supp <- function(dt, nm, tbl, filt) {
    write_source_data(dt, paste0("figureS_local_control_audit_unbounded", arm, "_", nm),
                      runs_used, tbl, SCRIPT, filt, data_dir)
}
TBL <- sprintf("results/combined/local-genetic-control-%s-{region}-vmrs.tsv", cohort)
## Panel b spans two modules: the decile comes from module 02's score, the
## outcome from module 03's end-to-end OOF R-squared. Both tables and both sets
## of run IDs are named, which is what puts module 03 into the run manifest's
## `upstream_runs` -- 03_close_figure_run.R unions `source_run_id` over the
## source-data tables, so a run that no panel names is a run the manifest does
## not record (AGENTS.md 7.11, 9).
TBL_PRED <- sprintf(paste("02_local_genetic_variance results/combined/local-genetic-control-%s-{region}-vmrs.tsv",
                          "+ 03_local_snp_prediction results/combined/oof-prediction-%s-{region}-vmrs.tsv"),
                    cohort, cohort)
sd(conc, "panel_a", TBL, FILTER)
sd(dec_sum, "panel_b", TBL_PRED,
   paste(FILTER, "; deciles of local_snp_contribution_score; outcome is",
         "r2_pred_oof joined on vmr_id; negative values retained (AGENTS.md 7.3)"),
   runs = runs_used_pred)
sd(pairs_dt, "panel_c", TBL,
   paste(FILTER, "; loci matched across regions by widest genomic overlap"))
sd(ctx_sum, "panel_d",
   paste(TBL, "+ 01_vmr_catalog qc/distance_to_nearest_gene.tsv"),
   paste(FILTER, "; joined on vmr_id; descriptive proportions, no enrichment test"))
sd_denom(excl_long, "panel_a", TBL, "all rows; eligibility as recorded upstream")
sd_denom(reasons, "panel_a_reasons", TBL, "local_genetic_control_eligible == FALSE")
sd_supp(elig[, .(n = .N, median = median(pve_cis_joint_unbounded),
            q25 = quantile(pve_cis_joint_unbounded, .25),
            q75 = quantile(pve_cis_joint_unbounded, .75)), by = region],
   "panel_a", TBL,
   paste(FILTER, "; AUDIT ONLY, not interpretable as absolute PVE"))

message("[done] Figure 2 written to ", fig_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
