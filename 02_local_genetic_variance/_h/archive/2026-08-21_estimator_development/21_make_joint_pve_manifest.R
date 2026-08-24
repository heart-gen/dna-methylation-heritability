#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))
source(file.path(dirname(script_path), "bslmm_pilot_functions.R"))
source(file.path(dirname(script_path), "joint_pve_functions.R"))

cli <- parse_cli(list(
    config = file.path(dirname(script_path), "..", "config",
                       "joint-pve-20260820.tsv"),
    mode = "validation",
    development_manifest = "",
    output = ""
))
if (!nzchar(cli$output)) stop("--output is required")
settings <- read_joint_settings(cli$config)
mode <- tolower(cli$mode)
if (!mode %in% c("development", "validation")) {
    stop("--mode must be development or validation")
}

required <- c(
    "sample_sizes", "num_snps", "ld_rho", "architectures", "h2_values",
    "validation_replicates_per_stratum", "base_seed", "validation_seed_offset",
    "outer_folds", "outer_repeats", "inner_folds", "alpha_grid",
    "lambda_rule", "max_features", "gemma_bin", "bslmm_mode",
    "bslmm_burn_in", "bslmm_sampling", "bslmm_rpace"
)
missing <- setdiff(required, names(settings))
if (length(missing)) stop("Missing settings: ", paste(missing, collapse = ", "))

if (identical(mode, "development")) {
    if (!nzchar(cli$development_manifest) ||
        !file.exists(cli$development_manifest)) {
        stop("Development mode requires --development_manifest")
    }
    manifest <- read_tsv(cli$development_manifest)
    expected_replicates <- 5L
    if (!identical(sort(unique(as.integer(manifest$replicate))),
                   seq_len(expected_replicates))) {
        stop("Development source must contain replicates 1--5")
    }
    manifest$split <- "development"
    manifest$feature_mode <- "development_augment"
} else {
    architectures <- split_character(settings$architectures)
    h2_values <- split_numeric(settings$h2_values)
    design_cells <- expand.grid(
        architecture = architectures,
        true_h2 = h2_values,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    strata <- expand.grid(
        n = as.integer(split_numeric(settings$sample_sizes)),
        num_snps = as.integer(split_numeric(settings$num_snps)),
        ld_rho = split_numeric(settings$ld_rho),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    reps <- as_int(settings$validation_replicates_per_stratum,
                   "validation_replicates_per_stratum")
    if (reps %% nrow(design_cells) != 0L) {
        stop("Validation replicates must balance architecture-by-PVE cells")
    }
    base_seed <- as_int(settings$base_seed, "base_seed") +
        as_int(settings$validation_seed_offset, "validation_seed_offset")
    records <- vector("list", nrow(strata))
    for (i in seq_len(nrow(strata))) {
        set.seed(base_seed + i)
        cell <- rep(seq_len(nrow(design_cells)), length.out = reps)
        cell <- sample(cell, length(cell), replace = FALSE)
        occurrence <- ave(cell, cell, FUN = seq_along)
        records[[i]] <- data.frame(
            split = "validation",
            feature_mode = "full",
            n = strata$n[[i]],
            num_snps = strata$num_snps[[i]],
            ld_rho = strata$ld_rho[[i]],
            architecture = design_cells$architecture[cell],
            true_h2 = design_cells$true_h2[cell],
            replicate = occurrence,
            stringsAsFactors = FALSE
        )
    }
    manifest <- do.call(rbind, records)
    manifest$scenario_id <- seq_len(nrow(manifest))
    manifest$seed <- base_seed + manifest$scenario_id * 1009L
}

manifest$outer_folds <- as_int(settings$outer_folds, "outer_folds")
manifest$outer_repeats <- as_int(settings$outer_repeats, "outer_repeats")
manifest$inner_folds <- as_int(settings$inner_folds, "inner_folds")
manifest$alpha_grid <- settings$alpha_grid
manifest$lambda_rule <- settings$lambda_rule
manifest$max_features <- as_int(settings$max_features, "max_features")
manifest$gemma_bin <- settings$gemma_bin
manifest$bslmm_mode <- as_int(settings$bslmm_mode, "bslmm_mode")
manifest$bslmm_burn_in <- as_int(settings$bslmm_burn_in, "bslmm_burn_in")
manifest$bslmm_sampling <- as_int(settings$bslmm_sampling, "bslmm_sampling")
manifest$bslmm_rpace <- as_int(settings$bslmm_rpace, "bslmm_rpace")
manifest <- manifest[, c(
    "scenario_id", "split", "feature_mode", "seed", "n", "num_snps",
    "ld_rho", "architecture", "true_h2", "replicate", "outer_folds",
    "outer_repeats", "inner_folds", "alpha_grid", "lambda_rule",
    "max_features", "gemma_bin", "bslmm_mode", "bslmm_burn_in",
    "bslmm_sampling", "bslmm_rpace"
)]
if (anyDuplicated(manifest$scenario_id)) stop("Duplicate scenario_id")
if (anyDuplicated(manifest$seed)) stop("Duplicate simulation seed")
write_tsv(manifest, cli$output)
cat("Wrote", nrow(manifest), mode, "joint-PVE scenarios to",
    normalizePath(cli$output), "\n")
