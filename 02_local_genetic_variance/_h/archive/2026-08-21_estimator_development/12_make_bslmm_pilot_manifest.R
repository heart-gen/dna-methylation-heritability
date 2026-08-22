#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))
source(file.path(dirname(script_path), "bslmm_pilot_functions.R"))

cli <- parse_cli(list(
    config = file.path(dirname(script_path), "..", "config", "bslmm-pilot.tsv"),
    output = file.path(dirname(script_path), "..", "_m", "config",
                       "bslmm-pilot-scenarios.tsv")
))
settings <- read_pilot_settings(cli$config)
required <- c(
    "sample_sizes", "num_snps", "ld_rho", "architectures", "h2_values",
    "replicates_per_cell", "base_seed", "outer_folds", "outer_repeats",
    "inner_folds", "alpha_grid", "lambda_rule", "max_features", "raw_metric",
    "max_design_distance", "bslmm_burn_in", "bslmm_sampling", "bslmm_rpace",
    "bslmm_mode", "gemma_bin"
)
missing <- setdiff(required, names(settings))
if (length(missing)) {
    stop("Missing pilot settings: ", paste(missing, collapse = ", "))
}

architectures <- trimws(strsplit(settings$architectures, ",", fixed = TRUE)[[1L]])
grid <- expand.grid(
    n = as.integer(split_numeric(settings$sample_sizes)),
    num_snps = as.integer(split_numeric(settings$num_snps)),
    ld_rho = split_numeric(settings$ld_rho),
    architecture = factor(architectures, levels = architectures),
    true_h2 = split_numeric(settings$h2_values),
    replicate = seq_len(as_int(settings$replicates_per_cell, "replicates_per_cell")),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
grid$architecture <- as.character(grid$architecture)
base_seed <- as_int(settings$base_seed, "base_seed")
## Keep expand.grid order (config order), not alphabetical architecture order.
rownames(grid) <- NULL
grid$scenario_id <- seq_len(nrow(grid))
grid$seed <- base_seed + grid$scenario_id * 1009L
grid$outer_folds <- as_int(settings$outer_folds, "outer_folds")
grid$outer_repeats <- as_int(settings$outer_repeats, "outer_repeats")
grid$inner_folds <- as_int(settings$inner_folds, "inner_folds")
grid$alpha_grid <- settings$alpha_grid
grid$lambda_rule <- settings$lambda_rule
grid$max_features <- as_int(settings$max_features, "max_features")
grid$raw_metric <- settings$raw_metric
grid$max_design_distance <- as_num(settings$max_design_distance, "max_design_distance")
grid$bslmm_mode <- as_int(settings$bslmm_mode, "bslmm_mode")
grid$bslmm_burn_in <- as_int(settings$bslmm_burn_in, "bslmm_burn_in")
grid$bslmm_sampling <- as_int(settings$bslmm_sampling, "bslmm_sampling")
grid$bslmm_rpace <- as_int(settings$bslmm_rpace, "bslmm_rpace")
grid$gemma_bin <- settings$gemma_bin
grid <- grid[, c(
    "scenario_id", "seed", "n", "num_snps", "ld_rho", "architecture",
    "true_h2", "replicate", "outer_folds", "outer_repeats", "inner_folds",
    "alpha_grid", "lambda_rule", "max_features", "raw_metric",
    "max_design_distance", "bslmm_mode", "bslmm_burn_in", "bslmm_sampling",
    "bslmm_rpace", "gemma_bin"
)]
write_tsv(grid, cli$output)
cat("Wrote", nrow(grid), "paired pilot scenarios to", normalizePath(cli$output), "\n")
