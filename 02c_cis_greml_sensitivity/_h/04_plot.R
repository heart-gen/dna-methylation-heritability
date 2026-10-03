#!/usr/bin/env Rscript
#### 02c -- run review figures ####
## Usage: Rscript _h/04_plot.R --run-id <id>

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02c_cis_greml_sensitivity", "_h")),
                 "cgs_functions.R"))

opts <- parse_v2_args(require = "run_id")
run_dir <- run_dir_for(opts$run_id)
run <- list(run_id = opts$run_id, dir = run_dir)
cfg <- load_cgs_config(run_dir)
man <- read_run_manifest(run_dir)
smoke <- identical(man[["smoke_run"]], "TRUE")
suppressPackageStartupMessages(library(ggplot2))
vmr <- fread(file.path(run_dir, "summary", "vmr-cis-greml.tsv"))
pf <- fread(file.path(run_dir, "summary", "decile-profile.tsv"))
sl <- fread(file.path(run_dir, "summary", "convergence-by-decile.tsv"))
save_both <- function(p, stem, w, h) {
    f <- file.path(run_dir, "figures", stem)
    grDevices::cairo_pdf(paste0(f, ".pdf"), width = w, height = h, bg = "white"); print(p); grDevices::dev.off()
    grDevices::png(paste0(f, ".png"), width = w, height = h, units = "in", res = 300,
                   bg = "white", type = "cairo"); print(p); grDevices::dev.off()
}
theme_set(theme_bw(base_size = 9))
d <- vmr[module02_eligible == TRUE & converged == TRUE]
p1 <- ggplot(d, aes(local_snp_contribution_score, h2_hat)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_point(size = 0.25, alpha = 0.2) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.5) +
    facet_wrap(~ reml_mode) +
    labs(x = "Module 02 local SNP contribution score (percentile)",
         y = "cis-GREML estimate (not reportable per VMR)",
         title = sprintf("%s, %s", man[["cohort"]], man[["region"]]))
save_both(p1, "cis_greml_vs_score", 8, 4)
p2 <- ggplot(pf, aes(score_decile, estimate, colour = reml_mode)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_pointrange(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(0.4),
                    size = 0.2) +
    scale_x_continuous(breaks = 1:10) +
    labs(x = "Module 02 score decile", y = "Mean cis-GREML estimate (jackknife 95% CI)")
p3 <- ggplot(sl, aes(score_decile, convergence_rate, colour = reml_mode)) +
    geom_line() + geom_point(size = 1) + scale_x_continuous(breaks = 1:10) +
    labs(x = "Module 02 score decile", y = "Converged fraction")
save_both(p2, "decile_profile", 6, 3.5)
save_both(p3, "convergence_by_decile", 6, 3)
