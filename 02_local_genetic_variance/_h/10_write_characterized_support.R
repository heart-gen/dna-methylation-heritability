#!/usr/bin/env Rscript

## Regenerate the frozen joint model's characterized feature support, PER CELL.
##
## Stage 03 previously bounded p_eff only by its mathematical range [1, n],
## because p_eff is a measured feature rather than a gridded design factor and
## the AR(1) training endpoints were accidental order statistics. That gate was
## blind: it passed all 11,239 eligible loci of lgv-AA-caudate-20260822 while
## 18.75% of them sat below the AR(1) minimum p_eff of 24.34 and 56.4%
## produced an unbounded estimate below the global minimum of every training
## simulation.
##
## The observed-regime grids characterise the same frozen model on real
## cis-window genotypes, one grid per cell x region. Taking the union across
## ALL cells -- which this script did until 2026-09-13 -- let one donor group
## widen another's eligibility domain: admitting the EA cells dropped the
## p_eff floor from 2.058 to 1.344, so an AA locus with p_eff = 1.5 became
## eligible on the strength of EA evidence alone. LD structure, MAF spectrum
## and effective rank differ between donor groups, so that borrowing has no
## biological justification.
##
## The support is therefore now computed SEPARATELY FOR EACH CELL, as the
## union of the AR(1) training grid (synthetic, ancestry-free, and so shared
## legitimately by every cell) and only that cell's own observed-regime grids.
## Grids are grouped by the `cohort` field of their own run manifest, which
## carries the cell token. Pooling across regions WITHIN a cell is retained:
## the concern is ancestry-driven LD, not tissue, and each region's production
## run needs its own donor count present in allowed_n.
##
## This script changes no estimator.  --regime-run-id takes a comma-separated
## list of run IDs.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
h_dir <- dirname(script_path)
source(file.path(h_dir, "00_functions.R"))

repo_root <- normalizePath(file.path(h_dir, "..", ".."))
module_root <- file.path(repo_root, "02_local_genetic_variance")
cli <- parse_cli(list(
    development_run_id = "lgv-joint-pve-train-20260820",
    regime_run_id = paste(c(
        ## arm cells (support the six accepted 20260823 arm production runs)
        "lgv-observed-regime-20260822",
        "lgv-observed-regime-AA-dlpfc-20260822",
        "lgv-observed-regime-AA-hippocampus-20260822",
        "lgv-observed-regime-all_individuals-caudate-20260822",
        "lgv-observed-regime-all_individuals-dlpfc-20260822",
        "lgv-observed-regime-all_individuals-hippocampus-20260822",
        ## donor-group estimation cells. The 20260910 (non-b) EA grids are
        ## superseded: they carried 2,388 computational failures from a
        ## mid-run working-tree mutation and must never characterise support.
        "lgv-observed-regime-all_individuals.AA-caudate-20260911",
        "lgv-observed-regime-all_individuals.AA-dlpfc-20260911",
        "lgv-observed-regime-all_individuals.AA-hippocampus-20260911",
        "lgv-observed-regime-all_individuals.EA-caudate-20260910b",
        "lgv-observed-regime-all_individuals.EA-dlpfc-20260910b",
        "lgv-observed-regime-all_individuals.EA-hippocampus-20260910b",
        ## tier-3 donor-count sensitivity cells (Module 08's caudate
        ## downsampling). A subsample cell needs its OWN grid: its donor count
        ## differs from the AA arm's, so `allowed_n` and the p_eff floor both
        ## move, and the arm's support would characterise a regime of 153
        ## donors that these cells never occupy.
        "lgv-observed-regime-AA.n118r1-caudate-20260918",
        "lgv-observed-regime-AA.n118r2-caudate-20260918",
        "lgv-observed-regime-AA.n118r3-caudate-20260918"
    ), collapse = ","),
    output = file.path(module_root, "config",
                       "joint-pve-characterized-support.tsv")
))
runs_root <- file.path(module_root, "_m", "runs")
development_path <- file.path(runs_root, cli$development_run_id, "combined",
                              "development-features.tsv")
regime_run_ids <- trimws(strsplit(cli$regime_run_id, ",", fixed = TRUE)[[1L]])
regime_run_ids <- regime_run_ids[nzchar(regime_run_ids)]
if (!length(regime_run_ids)) stop("No observed-regime run IDs supplied")
regime_paths <- file.path(runs_root, regime_run_ids, "results", "combined",
                          "observed-regime-estimates.tsv")
names(regime_paths) <- regime_run_ids
manifest_paths <- file.path(runs_root, regime_run_ids, "manifest.tsv")
names(manifest_paths) <- regime_run_ids
for (path in c(development_path, regime_paths, manifest_paths)) {
    if (!file.exists(path)) stop("Missing characterization input: ", path)
}
development <- read_tsv(development_path)
development <- development[development$feature_complete %in% TRUE, , drop = FALSE]
if (!nrow(development)) stop("Development grid is empty")

## Each grid declares its own cell in its run manifest; never infer it from the
## run ID, which is a naming convention rather than a recorded fact.
regime_cell <- vapply(regime_run_ids, function(run_id) {
    man <- read_tsv(manifest_paths[[run_id]])
    value <- man$value[man$field == "cohort"]
    if (length(value) != 1L) {
        stop("Regime run manifest lacks a unique cohort field: ", run_id)
    }
    as.character(value[[1L]])
}, character(1L))

regime_parts <- lapply(regime_run_ids, function(run_id) {
    part <- read_tsv(regime_paths[[run_id]])
    part <- part[part$feature_complete %in% TRUE, , drop = FALSE]
    if (!nrow(part)) stop("Observed-regime grid is empty: ", run_id)
    part <- part[, c("num_snps", "p_eff", "ld_metric", "n"), drop = FALSE]
    part$cell <- regime_cell[[run_id]]
    part
})
names(regime_parts) <- regime_run_ids

sha256 <- function(path) {
    tolower(sub(" .*$", "", system2("sha256sum", normalizePath(path),
                                    stdout = TRUE)[[1L]]))
}
development_sha <- sha256(development_path)
fields <- c("num_snps", "p_eff", "ld_metric")
cells <- sort(unique(regime_cell))

support <- do.call(rbind, lapply(cells, function(cell) {
    ids <- regime_run_ids[regime_cell == cell]
    regime <- do.call(rbind, regime_parts[ids])
    if (!nrow(regime)) stop("No regime rows for cell: ", cell)
    rows <- do.call(rbind, lapply(fields, function(field) {
        dev <- as.numeric(development[[field]])
        reg <- as.numeric(regime[[field]])
        data.frame(
            cell = cell,
            feature = field,
            development_min = min(dev, na.rm = TRUE),
            development_max = max(dev, na.rm = TRUE),
            regime_min = min(reg, na.rm = TRUE),
            regime_max = max(reg, na.rm = TRUE),
            support_min = min(c(dev, reg), na.rm = TRUE),
            support_max = max(c(dev, reg), na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    }))
    rows$development_run_id <- cli$development_run_id
    rows$regime_run_id <- paste(ids, collapse = ",")
    rows$development_sha256 <- development_sha
    rows$regime_sha256 <- paste(vapply(regime_paths[ids], sha256,
                                       character(1L)), collapse = ",")
    rows$allowed_n <- paste(sort(unique(c(as.integer(development$n),
                                          as.integer(regime$n)))),
                            collapse = ",")
    rows
}))
write_tsv(support, cli$output)
cat("Wrote per-cell characterized support:", length(cells), "cells x",
    length(fields), "features to", normalizePath(cli$output), "\n\n")
print(support[, c("cell", "feature", "support_min", "support_max")],
      row.names = FALSE)
cat("\nallowed_n by cell:\n")
print(unique(support[, c("cell", "allowed_n")]), row.names = FALSE)
