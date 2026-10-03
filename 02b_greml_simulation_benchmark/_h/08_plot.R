#!/usr/bin/env Rscript
#### 02b -- run figures ####
##
## Review figures for the run itself; the manuscript supplement is built by
## 11_integrated_manuscript_outputs from _m/combined/. PDF through cairo_pdf
## (base pdf() drops non-ASCII glyphs; see 11/_h/00_figure_theme.R).
##
## Usage: Rscript _h/08_plot.R --run-id <id>

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02b_greml_simulation_benchmark", "_h")),
                 "greml_functions.R"))
suppressPackageStartupMessages(library(ggplot2))

opts <- parse_v2_args(require = "run_id")
run_dir <- run_dir_for(opts$run_id)
cfg <- load_greml_config(run_dir)
man <- read_run_manifest(run_dir)
arm1 <- identical(man[["arm"]], cfg$arm1$id)
est <- fread(file.path(run_dir, "summary", "estimates.tsv"))
met <- fread(file.path(run_dir, "summary", "recovery-metrics.tsv"))

save_both <- function(p, stem, w, h) {
    f <- file.path(run_dir, "figures", stem)
    grDevices::cairo_pdf(paste0(f, ".pdf"), width = w, height = h, bg = "white")
    print(p); grDevices::dev.off()
    grDevices::png(paste0(f, ".png"), width = w, height = h, units = "in", res = 300,
                   bg = "white", type = "cairo")
    print(p); grDevices::dev.off()
}
theme_set(theme_bw(base_size = 9))

if (arm1) {
    est[, n_label := factor(paste0("n = ", N), levels = paste0("n = ", sort(unique(N))))]
    p1 <- ggplot(est, aes(truth_h2, h2_hat)) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
        geom_point(size = 0.3, alpha = 0.3) +
        geom_smooth(method = "loess", se = FALSE, linewidth = 0.5, formula = y ~ x) +
        facet_grid(reml_mode ~ n_label) +
        labs(x = "Realized simulated h2", y = "GREML-LDMS estimate",
             title = "Arm 1: AR(1) simulated genotypes (out of regime)")
    save_both(p1, "arm1_estimate_vs_truth", 10, 4.5)
    m <- melt(met, id.vars = c("cell", "reml_mode"),
              measure.vars = c("bias", "rmse", "coverage95", "estimation_failure_rate"))
    p2 <- ggplot(m, aes(cell, value, colour = reml_mode)) +
        geom_line() + geom_point(size = 1) + scale_x_log10() +
        facet_wrap(~variable, scales = "free_y", nrow = 1) +
        labs(x = "Simulated sample size (log scale)", y = NULL)
    save_both(p2, "arm1_metrics_by_n", 10, 3)
} else {
    p1 <- ggplot(est, aes(factor(h2_nominal), h2_hat)) +
        geom_hline(yintercept = c(0, 1), colour = "grey80") +
        geom_boxplot(outlier.size = 0.2, linewidth = 0.3) +
        stat_summary(aes(y = truth_h2, group = 1), fun = mean, geom = "line",
                     colour = "red", linewidth = 0.4) +
        facet_grid(reml_mode ~ architecture) +
        labs(x = "Simulated h2", y = "Cis-window REML estimate",
             title = sprintf("Arm 2: real AA cis-window genotypes, %s, n = %s",
                             man[["region"]], man[["design_n"]]))
    save_both(p1, "arm2_estimate_by_h2", 8, 5)
    m <- melt(met, id.vars = c("reml_mode", "architecture", "h2_nominal"),
              measure.vars = c("bias", "rmse", "coverage95", "estimation_failure_rate"))
    p2 <- ggplot(m, aes(h2_nominal, value, colour = architecture)) +
        geom_line() + geom_point(size = 1) +
        facet_grid(reml_mode ~ variable, scales = "free_y") +
        labs(x = "Simulated h2", y = NULL)
    save_both(p2, "arm2_metrics_by_h2", 10, 4)
}
