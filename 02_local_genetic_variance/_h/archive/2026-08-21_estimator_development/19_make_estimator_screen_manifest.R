#!/usr/bin/env Rscript

## Estimator-screening manifest. Every arm receives the same design cells with
## the same seed, so simulate_locus produces byte-identical draws and the arms
## are exactly paired. Only outer_repeats, alpha_grid, and lambda_rule vary.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    config = file.path(dirname(script_path), "..", "config",
                       "estimator-screen-20260818.tsv"),
    arms = file.path(dirname(script_path), "..", "config",
                     "estimator-screen-arms.tsv"),
    output = ""
))
if (!nzchar(cli$output)) stop("--output is required")

config_table <- read_tsv(cli$config)
settings <- stats::setNames(as.list(config_table$value), config_table$field)
arms <- read_tsv(cli$arms)
required_arm <- c("arm_id", "outer_repeats", "alpha_grid", "lambda_rule")
missing_arm <- setdiff(required_arm, names(arms))
if (length(missing_arm)) {
    stop("Arm table missing columns: ", paste(missing_arm, collapse = ", "))
}
if (anyDuplicated(arms$arm_id)) stop("Arm identifiers must be unique")

sample_sizes <- as.integer(split_numeric(settings$sample_sizes))
num_snps <- as.integer(split_numeric(settings$num_snps))
ld_rho <- split_numeric(settings$ld_rho)
h2_values <- split_numeric(settings$h2_values)
architectures <- trimws(strsplit(settings$architectures, ",", fixed = TRUE)[[1L]])
replicates <- as_int(settings$replicates_per_cell, "replicates_per_cell")
base_seed <- as_int(settings$base_seed, "base_seed") +
    as_int(settings$seed_offset, "seed_offset")

cells <- expand.grid(
    replicate = seq_len(replicates),
    true_h2 = h2_values,
    architecture = architectures,
    ld_rho = ld_rho,
    num_snps = num_snps,
    n = sample_sizes,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
cells <- cells[order(cells$n, cells$num_snps, cells$ld_rho,
                     cells$architecture, cells$true_h2, cells$replicate), ]
## The seed identifies the design cell only, never the arm, so matched cells
## across arms simulate the same locus.
cells$cell_id <- seq_len(nrow(cells))
cells$seed <- base_seed + cells$cell_id * 1009L

records <- lapply(seq_len(nrow(arms)), function(i) {
    out <- cells
    out$split <- arms$arm_id[[i]]
    out$arm_id <- arms$arm_id[[i]]
    out$outer_repeats <- as.integer(arms$outer_repeats[[i]])
    out$alpha_grid <- arms$alpha_grid[[i]]
    out$lambda_rule <- arms$lambda_rule[[i]]
    out
})
manifest <- do.call(rbind, records)
manifest$scenario_id <- seq_len(nrow(manifest))
manifest$outer_folds <- as_int(settings$outer_folds, "outer_folds")
manifest$inner_folds <- as_int(settings$inner_folds, "inner_folds")
manifest$max_features <- as_int(settings$max_features, "max_features")
manifest$raw_metric <- settings$raw_metric
manifest$null_alpha <- 0.05
manifest$max_design_distance <- 1.25

manifest <- manifest[, c(
    "scenario_id", "split", "arm_id", "cell_id", "seed", "n", "num_snps",
    "ld_rho", "architecture", "true_h2", "replicate", "outer_folds",
    "outer_repeats", "inner_folds", "alpha_grid", "lambda_rule",
    "max_features", "raw_metric", "null_alpha", "max_design_distance"
)]
write_tsv(manifest, cli$output)

cat("Wrote", nrow(manifest), "scenarios to", normalizePath(cli$output), "\n")
cat("Arms:", paste(arms$arm_id, collapse = ", "), "\n")
cat("Design cells per arm:", nrow(cells), "\n")
stopifnot(
    "each arm must cover every design cell" =
        all(table(manifest$arm_id) == nrow(cells)),
    "matched cells must share a seed across arms" =
        all(tapply(manifest$seed, manifest$cell_id,
                   function(x) length(unique(x))) == 1L)
)
cat("Pairing check passed: matched cells share one seed across all arms\n")
