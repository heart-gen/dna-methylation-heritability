#!/usr/bin/env Rscript

## Regenerate the frozen joint model's characterized feature support.
##
## Stage 03 previously bounded p_eff only by its mathematical range [1, n],
## because p_eff is a measured feature rather than a gridded design factor and
## the AR(1) training endpoints were accidental order statistics. That gate was
## blind: it passed all 11,239 eligible loci of lgv-AA-caudate-20260822 while
## 18.75% of them sat below the AR(1) minimum p_eff of 24.34 and 56.4%
## produced an unbounded estimate below the global minimum of every training
## simulation.
##
## The 2026-08-22 observed-regime grid characterises the same frozen model on
## real cis-window genotypes, so the defensible support is the union of the two
## grids. This script writes that union; it changes no estimator.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
h_dir <- dirname(script_path)
source(file.path(h_dir, "00_functions.R"))

repo_root <- normalizePath(file.path(h_dir, "..", ".."))
module_root <- file.path(repo_root, "02_local_genetic_variance")
cli <- parse_cli(list(
    development_run_id = "lgv-joint-pve-train-20260820",
    regime_run_id = "lgv-observed-regime-20260822",
    output = file.path(module_root, "config",
                       "joint-pve-characterized-support.tsv")
))
runs_root <- file.path(module_root, "_m", "runs")
development_path <- file.path(runs_root, cli$development_run_id, "combined",
                              "development-features.tsv")
regime_path <- file.path(runs_root, cli$regime_run_id, "results", "combined",
                         "observed-regime-estimates.tsv")
for (path in c(development_path, regime_path)) {
    if (!file.exists(path)) stop("Missing characterization input: ", path)
}
development <- read_tsv(development_path)
development <- development[development$feature_complete %in% TRUE, , drop = FALSE]
regime <- read_tsv(regime_path)
regime <- regime[regime$feature_complete %in% TRUE, , drop = FALSE]
if (!nrow(development) || !nrow(regime)) stop("A characterization grid is empty")

sha256 <- function(path) {
    tolower(sub(" .*$", "", system2("sha256sum", normalizePath(path),
                                    stdout = TRUE)[[1L]]))
}
fields <- c("num_snps", "p_eff", "ld_metric")
records <- lapply(fields, function(field) {
    dev <- as.numeric(development[[field]])
    reg <- as.numeric(regime[[field]])
    data.frame(
        feature = field,
        development_min = min(dev, na.rm = TRUE),
        development_max = max(dev, na.rm = TRUE),
        regime_min = min(reg, na.rm = TRUE),
        regime_max = max(reg, na.rm = TRUE),
        support_min = min(c(dev, reg), na.rm = TRUE),
        support_max = max(c(dev, reg), na.rm = TRUE),
        stringsAsFactors = FALSE
    )
})
support <- do.call(rbind, records)
support$development_run_id <- cli$development_run_id
support$regime_run_id <- cli$regime_run_id
support$development_sha256 <- sha256(development_path)
support$regime_sha256 <- sha256(regime_path)
support$allowed_n <- paste(sort(unique(c(as.integer(development$n),
                                         as.integer(regime$n)))),
                           collapse = ",")
write_tsv(support, cli$output)
cat("Wrote characterized support for", nrow(support), "features to",
    normalizePath(cli$output), "\n")
print(support[, c("feature", "support_min", "support_max")])
