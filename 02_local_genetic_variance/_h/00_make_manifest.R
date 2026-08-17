#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    config = file.path(dirname(script_path), "..", "config", "analysis.tsv"),
    output = file.path(dirname(script_path), "..", "_m", "config", "scenarios.tsv"),
    split = "both",
    seed_offset = "0"
))

config_table <- read_tsv(cli$config)
if (!all(c("setting", "value") %in% names(config_table))) {
    stop("Configuration must contain setting and value columns")
}
settings <- stats::setNames(as.list(config_table$value), config_table$setting)
required <- c(
    "sample_sizes", "num_snps", "ld_rho", "architectures", "h2_values",
    "calibration_replicates_per_stratum", "evaluation_replicates_per_stratum",
    "outer_folds", "outer_repeats", "inner_folds", "alpha_grid",
    "lambda_rule", "max_features", "base_seed", "raw_metric",
    "max_design_distance"
)
missing <- setdiff(required, names(settings))
if (length(missing)) stop("Missing configuration settings: ", paste(missing, collapse = ", "))
if (!"null_alpha" %in% names(settings)) {
    warning("Legacy configuration lacks null_alpha; using 0.05")
    settings$null_alpha <- "0.05"
}

sample_sizes <- as.integer(split_numeric(settings$sample_sizes))
num_snps <- as.integer(split_numeric(settings$num_snps))
ld_rho <- split_numeric(settings$ld_rho)
h2_values <- split_numeric(settings$h2_values)
architectures <- strsplit(settings$architectures, ",", fixed = TRUE)[[1L]]
architectures <- trimws(architectures)
base_seed <- as_int(settings$base_seed, "base_seed") +
    as_int(cli$seed_offset, "seed_offset")

strata <- expand.grid(
    n = sample_sizes,
    num_snps = num_snps,
    ld_rho = ld_rho,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
design_cells <- expand.grid(
    architecture = architectures,
    true_h2 = h2_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)

make_split <- function(split_name, replicates_per_stratum, seed_offset) {
    if (replicates_per_stratum < nrow(design_cells) ||
        replicates_per_stratum %% nrow(design_cells) != 0L) {
        stop(
            split_name, " replicates_per_stratum must be a positive multiple of ",
            nrow(design_cells),
            " so every architecture-by-h2 cell is exactly balanced"
        )
    }
    records <- vector("list", nrow(strata))
    for (i in seq_len(nrow(strata))) {
        set.seed(base_seed + seed_offset + i)
        cell_index <- rep(seq_len(nrow(design_cells)), length.out = replicates_per_stratum)
        cell_index <- sample(cell_index, length(cell_index), replace = FALSE)
        cell_occurrence <- ave(cell_index, cell_index, FUN = seq_along)
        records[[i]] <- data.frame(
            split = split_name,
            n = strata$n[[i]],
            num_snps = strata$num_snps[[i]],
            ld_rho = strata$ld_rho[[i]],
            architecture = design_cells$architecture[cell_index],
            true_h2 = design_cells$true_h2[cell_index],
            replicate = cell_occurrence,
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, records)
}

calibration <- make_split(
    "calibration",
    as_int(settings$calibration_replicates_per_stratum,
           "calibration_replicates_per_stratum"),
    100000L
)
evaluation <- make_split(
    "evaluation",
    as_int(settings$evaluation_replicates_per_stratum,
           "evaluation_replicates_per_stratum"),
    200000L
)
manifest <- rbind(calibration, evaluation)
requested_split <- tolower(cli$split)
if (!requested_split %in% c("both", "calibration", "evaluation")) {
    stop("split must be both, calibration, or evaluation")
}
if (requested_split != "both") {
    manifest <- manifest[manifest$split == requested_split, , drop = FALSE]
}
manifest$scenario_id <- seq_len(nrow(manifest))
manifest$seed <- base_seed + manifest$scenario_id * 1009L
manifest$outer_folds <- as_int(settings$outer_folds, "outer_folds")
manifest$outer_repeats <- as_int(settings$outer_repeats, "outer_repeats")
manifest$inner_folds <- as_int(settings$inner_folds, "inner_folds")
manifest$alpha_grid <- settings$alpha_grid
manifest$lambda_rule <- settings$lambda_rule
manifest$max_features <- as_int(settings$max_features, "max_features")
manifest$raw_metric <- settings$raw_metric
manifest$null_alpha <- as_num(settings$null_alpha, "null_alpha")
manifest$max_design_distance <- as_num(settings$max_design_distance, "max_design_distance")
manifest <- manifest[, c(
    "scenario_id", "split", "seed", "n", "num_snps", "ld_rho",
    "architecture", "true_h2", "replicate", "outer_folds",
    "outer_repeats", "inner_folds", "alpha_grid", "lambda_rule",
    "max_features", "raw_metric", "null_alpha", "max_design_distance"
)]
write_tsv(manifest, cli$output)

cat("Wrote", nrow(manifest), "scenarios to", normalizePath(cli$output), "\n")
cat("Calibration:", sum(manifest$split == "calibration"), "\n")
cat("Evaluation:", sum(manifest$split == "evaluation"), "\n")
