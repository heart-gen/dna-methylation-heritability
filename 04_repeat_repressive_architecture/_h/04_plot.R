#!/usr/bin/env Rscript
#### 04_repeat_repressive_architecture -- figures for the cross-region result ####
##
## Runs once per cohort, after 03_apply_gates.R, and reads only that stage's
## pooled outputs. Every panel is drawn on the PRIMARY analysis set unless the
## panel is explicitly about a sensitivity, so a reader cannot mistake a
## sensitivity estimate for the headline one.
##
## The estimates are on the log-odds scale, per 1 SD of the within-cell local
## SNP contribution score. They are NOT PVE differences and the axis labels say
## so, because a figure travels further than its caption.
##
##   Rscript _h/04_plot.R --cohort AA --run-ids id1,id2,id3

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages({
    library(data.table); library(ggplot2)
})

MODULE <- "04_repeat_repressive_architecture"
opts <- parse_v2_args(require = c("cohort", "run_ids"))
run_ids <- trimws(strsplit(opts$run_ids, ",")[[1]])
out_dir <- file.path(repo_root(), MODULE, "_m", "runs", run_ids[1], "results")

res_f <- file.path(out_dir, "association-results-all-regions.tsv")
claims_f <- file.path(out_dir, "interpretation-claims.tsv")
for (f in c(res_f, claims_f)) {
    if (!file.exists(f)) stop("Missing ", f, "; run 03_apply_gates.R first.")
}
res <- fread(res_f)
claims <- fread(claims_f)

fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(p, name, w = 7, h = 4.5) {
    ggsave(file.path(fig_dir, paste0(name, ".png")), p,
           width = w, height = h, dpi = 300)
}

cfg <- load_config("repeat_annotations")
alpha <- config_get(cfg, "interpretation.alpha")

## ------------------------------------------------------- 1. forest, primary
## One row per region x outcome with a 95% Wald interval. The zero line is the
## null; a claim is permitted only where the required number of regions clear
## it in the SAME direction, so drawing all regions together is the point.
prim <- res[analysis_set == "primary" &
                predictor == config_get(cfg, "primary_model.predictor")]
if (nrow(prim) > 0) {
    prim[, `:=`(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)]
    p <- ggplot(prim, aes(x = estimate, y = region, colour = p < alpha)) +
        geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
        geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.18) +
        geom_point(size = 2) +
        facet_wrap(~ outcome, ncol = 1, scales = "free_x") +
        scale_colour_manual(values = c(`TRUE` = "#1b6ca8", `FALSE` = "grey55"),
                            name = sprintf("p < %s", alpha)) +
        labs(x = paste("log-odds per 1 SD of the local SNP contribution score",
                       "(NOT a PVE difference)"),
             y = NULL,
             title = paste0("Repressive-compartment association, ",
                            opts$cohort, ", primary analysis set")) +
        theme_bw(base_size = 10)
    save_fig(p, "forest-primary", h = 2 + 1.6 * uniqueN(prim$outcome))
}

## ------------------------------------- 2. sensitivity stability per outcome
## The gates require an estimate to survive every locked sensitivity, so the
## sensitivities belong in the figure set rather than in a supplement nobody
## opens. Colour marks the primary set so it stays distinguishable.
sens <- res[predictor == config_get(cfg, "primary_model.predictor")]
if (uniqueN(sens$analysis_set) > 1) {
    sens[, `:=`(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)]
    p <- ggplot(sens, aes(x = estimate, y = analysis_set,
                          colour = analysis_set == "primary")) +
        geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
        geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.18) +
        geom_point(size = 1.8) +
        facet_grid(outcome ~ region, scales = "free_x") +
        scale_colour_manual(values = c(`TRUE` = "#1b6ca8", `FALSE` = "grey55"),
                            guide = "none") +
        labs(x = "log-odds per 1 SD (NOT a PVE difference)", y = NULL,
             title = paste0("Sensitivity stability, ", opts$cohort)) +
        theme_bw(base_size = 9)
    save_fig(p, "sensitivity-stability",
             w = 3 + 2.4 * uniqueN(sens$region),
             h = 2 + 1.8 * uniqueN(sens$outcome))
}

## --------------------------------------------------- 3. what may be claimed
## The claims table is the module's actual conclusion, so it is rendered as a
## panel rather than left as a TSV. Regions that do not survive are drawn, not
## omitted -- a blank cell is evidence.
## Selected by intersection with the regions actually present, NOT by
## subtracting a hardcoded list of metadata columns: the subtraction form
## silently treats any new claims column (e.g. concentration_supported) as a
## region and plots it as one.
region_cols <- intersect(names(claims), unique(prim$region))
if (length(region_cols) > 0) {
    long <- melt(claims, id.vars = "outcome", measure.vars = region_cols,
                 variable.name = "region", value.name = "survives")
    long[, survives := as.logical(survives)]
    p <- ggplot(long, aes(x = region, y = outcome, fill = survives)) +
        geom_tile(colour = "white", linewidth = 1) +
        scale_fill_manual(values = c(`TRUE` = "#1b6ca8", `FALSE` = "grey85"),
                          na.value = "grey95", name = "survives gates") +
        labs(x = NULL, y = NULL,
             title = paste0("Regions surviving the locked sensitivities, ",
                            opts$cohort),
             caption = "Overlap only; not repeat activity, expression, or retrotransposition") +
        theme_minimal(base_size = 10)
    save_fig(p, "gate-survival", w = 5.5, h = 1.4 + 0.5 * uniqueN(long$outcome))
}

message("[04] figures written to ", fig_dir)
