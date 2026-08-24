#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    input = file.path(dirname(script_path), "..", "_m", "raw", "evaluation"),
    model = file.path(dirname(script_path), "..", "_m", "calibration", "elastic-net-calibration.rds"),
    analysis_config = file.path(dirname(script_path), "..", "config", "analysis.tsv"),
    output_dir = file.path(dirname(script_path), "..", "_m", "evaluation")
))
analysis_config <- read_tsv(cli$analysis_config)
analysis_settings <- stats::setNames(analysis_config$value, analysis_config$setting)
if (!"low_h2_max" %in% names(analysis_settings)) {
    stop("Analysis configuration must lock low_h2_max")
}
low_h2_max <- as_num(analysis_settings[["low_h2_max"]], "low_h2_max")
files <- list.files(cli$input, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)
if (!length(files)) stop("No evaluation result files found in ", cli$input)
evaluation <- do.call(rbind, lapply(files, read_tsv))
evaluation <- evaluation[evaluation$converged %in% TRUE, , drop = FALSE]
model <- readRDS(cli$model)
calibrated <- predict_calibration(model, evaluation)
evaluation <- cbind(evaluation, calibrated)
evaluation$error <- evaluation$h2_en_calibrated - evaluation$true_h2
evaluation$he_error <- evaluation$he_h2 - evaluation$true_h2
evaluation$covered <- evaluation$true_h2 >= evaluation$h2_calibration_lower &
    evaluation$true_h2 <= evaluation$h2_calibration_upper

summarize_group <- function(data) {
    available <- is.finite(data$h2_en_calibrated)
    mean_or_na <- function(x) if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_
    data.frame(
        simulations = nrow(data),
        estimable_simulations = sum(available),
        estimable_rate = mean(available),
        mean_true_h2 = mean_or_na(data$true_h2),
        mean_estimated_h2 = mean_or_na(data$h2_en_calibrated),
        bias = mean_or_na(data$error),
        absolute_mean_bias = abs(mean_or_na(data$error)),
        mean_absolute_error = mean_or_na(abs(data$error)),
        rmse = sqrt(mean_or_na(data$error^2)),
        calibration_interval_coverage = mean_or_na(data$covered),
        positive_signal_rate = mean_or_na(data$positive_signal),
        within_domain_rate = mean(data$calibration_status == "within_domain"),
        he_convergence_rate = mean(data$he_converged %in% TRUE),
        he_bias = mean_or_na(data$he_error),
        he_rmse = sqrt(mean_or_na(data$he_error^2)),
        he_out_of_bounds_rate = mean_or_na(data$he_h2 < 0 | data$he_h2 > 1),
        calibrated_out_of_bounds_rate = mean_or_na(
            data$h2_en_calibrated < 0 | data$h2_en_calibrated > 1
        ),
        stringsAsFactors = FALSE
    )
}

group_columns <- c("n", "num_snps", "ld_rho", "architecture", "true_h2")
group_key <- interaction(evaluation[group_columns], drop = TRUE, lex.order = TRUE)
by_design <- do.call(rbind, lapply(split(evaluation, group_key), function(data) {
    cbind(data[1L, group_columns, drop = FALSE], summarize_group(data))
}))
rownames(by_design) <- NULL
overall <- summarize_group(evaluation)
overall$spearman_truth_estimate <- suppressWarnings(stats::cor(
    evaluation$true_h2, evaluation$h2_en_calibrated,
    method = "spearman", use = "complete.obs"
))
overall$he_spearman_truth_estimate <- suppressWarnings(stats::cor(
    evaluation$true_h2, evaluation$he_h2,
    method = "spearman", use = "complete.obs"
))
null <- evaluation[evaluation$true_h2 == 0, , drop = FALSE]
overall$null_type1_error <- if (nrow(null)) mean(null$positive_signal) else NA_real_
h2_level_bias <- aggregate(error ~ true_h2, evaluation, mean)
low <- evaluation$true_h2 <= low_h2_max
low_h2_level_bias <- h2_level_bias[
    h2_level_bias$true_h2 <= low_h2_max, , drop = FALSE
]
if (!any(low) || !nrow(low_h2_level_bias)) {
    stop("No evaluation simulations fall at or below low_h2_max=", low_h2_max)
}
overall$null_mean_estimated_h2 <- if (nrow(null)) {
    mean(null$h2_en_calibrated)
} else {
    NA_real_
}
overall$low_h2_max <- low_h2_max
overall$low_h2_absolute_mean_bias <- abs(mean(evaluation$error[low]))
overall$max_absolute_low_h2_level_bias <- max(abs(low_h2_level_bias$error))
overall$max_absolute_h2_level_bias <- max(abs(h2_level_bias$error))

dir.create(cli$output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(evaluation, file.path(cli$output_dir, "calibrated-evaluation-estimates.tsv"))
write_tsv(by_design, file.path(cli$output_dir, "calibration-performance-by-design.tsv"))
write_tsv(h2_level_bias, file.path(cli$output_dir, "calibration-bias-by-h2.tsv"))
write_tsv(overall, file.path(cli$output_dir, "calibration-performance-overall.tsv"))
capture_session_info(file.path(cli$output_dir, "session-info.txt"))
cat("Wrote evaluation outputs to", normalizePath(cli$output_dir), "\n")
