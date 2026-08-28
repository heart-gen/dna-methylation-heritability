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
##   A  the rank is reproduced by independent estimators, and is NOT explained
##      by locus geometry (the reviewer's first objection)
##   B  held-out prediction accuracy rises monotonically across the rank
##      -- the secondary endpoint of AGENTS.md 11
##   C  the rank is concordant across brain regions in VMRs called in both
##   D  genic context across the rank -- descriptive proportions only, no
##      enrichment test, no threshold. Repeat and repressive-chromatin
##      enrichment is Figure 3 / Module 04, whose enrichment model is a locked
##      PI decision (AGENTS.md 12); nothing here anticipates it.
##   E  denominators and exclusions (AGENTS.md 7.9)
##
## Usage:
##   Rscript 02_figure2_local_control.R --cohort AA --run-id fig-all-20260826

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(V2_ROOT, "10_integrated_manuscript_outputs", "_h",
                 "00_figure_theme.R"))

suppressPackageStartupMessages({
    library(patchwork)
    library(scales)
    library(GenomicRanges)
})

opts <- parse_v2_args(require = c("cohort", "run_id"))
cohort <- opts$cohort
regions <- load_config("cohorts")$regions

module_root <- file.path(V2_ROOT, "10_integrated_manuscript_outputs")
run_dir  <- file.path(module_root, "_m", "runs", opts$run_id)
fig_dir  <- file.path(run_dir, "figures")
data_dir <- file.path(run_dir, "source_data")

LGV_RUN <- function(r) paste0("lgv-", cohort, "-", r, "-20260823")
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

SCRIPT <- "10_integrated_manuscript_outputs/_h/02_figure2_local_control.R"
runs_used <- vapply(regions, LGV_RUN, "")
FILTER <- "local_genetic_control_eligible == TRUE"

## ------------------------- A. what the rank agrees with, and what it does not
##
## Two contrasting groups on one axis. Independent estimators of local genetic
## control should track the rank; locus geometry should not, because a ranking
## driven by SNP count or LD would be an artifact rather than a signal.
CONCORD <- c(bslmm_pve = "BSLMM PVE",
             he_h2     = "Haseman-Elston",
             rho2_oof  = "Out-of-fold \u03c1\u00b2",
             r2_oof    = "Held-out R\u00b2")
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
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                   position = position_dodge(width = 0.62), linewidth = 0.42) +
    geom_point(size = 1.5, position = position_dodge(width = 0.62)) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = REGION_COLORS, name = NULL) +
    scale_x_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Spearman correlation with local SNP contribution rank", y = NULL) +
    BASE_THEME + NO_TITLES +
    theme(legend.position = "top", legend.margin = margin(0, 0, -4, 0),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          strip.text.y.right = element_text(angle = -90, face = "bold"),
          plot.margin = margin(5, 12, 5, 8))

## ------------------------------ B. held-out prediction across the rank
##
## The secondary endpoint. Deciles of the rank, not thresholds: no cut point is
## claimed, and the axis is prediction accuracy, never variance explained.
dec <- copy(elig)
dec[, decile := cut(local_snp_contribution_score, breaks = seq(0, 1, 0.1),
                    labels = 1:10, include.lowest = TRUE)]
dec_sum <- dec[!is.na(decile), .(
    n = .N,
    median = median(r2_oof, na.rm = TRUE),
    q25 = quantile(r2_oof, 0.25, na.rm = TRUE),
    q75 = quantile(r2_oof, 0.75, na.rm = TRUE)), by = .(region, decile)]

pB <- ggplot(dec_sum, aes(as.integer(decile), median, colour = region, fill = region)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.3) +
    scale_colour_manual(values = REGION_COLORS, guide = "none") +
    scale_fill_manual(values = REGION_COLORS, guide = "none") +
    scale_x_continuous(breaks = 1:10) +
    labs(x = "Decile of local SNP contribution rank",
         y = expression("Held-out"~R^2)) +
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
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
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
CTX_RUN <- function(r) paste0("vmrcatqc-", cohort, "-", r, "-20260826-a")
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

pD <- ggplot(excl_long, aes(region, n, fill = status)) +
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
figure <- (pA / (pB | pC) / pCtx / pD) +
    plot_layout(heights = c(1.25, 1.0, 0.85, 0.85)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 11))

save_figure(figure, paste0("figure2_local_genetic_control", arm),
            width = FIG_WIDTH_FULL, height = 10.2, fig_dir = fig_dir)

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
sd <- function(dt, nm, tbl, filt) {
    write_source_data(dt, paste0("figure2_local_genetic_control", arm, "_", nm),
                      runs_used, tbl, SCRIPT, filt, data_dir)
}
## The audit panel ships as its own supplemental figure, so its source data is
## named for that figure rather than for Figure 2.
sd_supp <- function(dt, nm, tbl, filt) {
    write_source_data(dt, paste0("figureS_local_control_audit_unbounded", arm, "_", nm),
                      runs_used, tbl, SCRIPT, filt, data_dir)
}
TBL <- sprintf("results/combined/local-genetic-control-%s-{region}-vmrs.tsv", cohort)
sd(conc, "panelA", TBL, FILTER)
sd(dec_sum, "panelB", TBL, paste(FILTER, "; deciles of local_snp_contribution_score"))
sd(pairs_dt, "panelC", TBL,
   paste(FILTER, "; loci matched across regions by widest genomic overlap"))
sd(ctx_sum, "panelD_genic_context",
   paste(TBL, "+ 01_vmr_catalog qc/distance_to_nearest_gene.tsv"),
   paste(FILTER, "; joined on vmr_id; descriptive proportions, no enrichment test"))
sd(excl_long, "panelE", TBL, "all rows; eligibility as recorded upstream")
sd(reasons, "panelE_reasons", TBL, "local_genetic_control_eligible == FALSE")
sd_supp(elig[, .(n = .N, median = median(pve_cis_joint_unbounded),
            q25 = quantile(pve_cis_joint_unbounded, .25),
            q75 = quantile(pve_cis_joint_unbounded, .75)), by = region],
   "distribution", TBL,
   paste(FILTER, "; AUDIT ONLY, not interpretable as absolute PVE"))

message("[done] Figure 2 written to ", fig_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
