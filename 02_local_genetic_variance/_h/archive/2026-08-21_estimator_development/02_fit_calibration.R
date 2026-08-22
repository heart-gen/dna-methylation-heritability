#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    input = file.path(dirname(script_path), "..", "_m", "raw", "calibration"),
    criteria = file.path(dirname(script_path), "..", "config", "acceptance-criteria.tsv"),
    analysis_config = file.path(dirname(script_path), "..", "config", "analysis.tsv"),
    output_model = file.path(dirname(script_path), "..", "_m", "calibration", "elastic-net-calibration.rds"),
    output_manifest = file.path(dirname(script_path), "..", "_m", "calibration", "calibration-manifest.tsv"),
    output_tuning = file.path(dirname(script_path), "..", "_m", "calibration", "hybrid-weight-tuning.tsv"),
    session_info = file.path(dirname(script_path), "..", "_m", "calibration", "session-info.txt")
))

files <- list.files(cli$input, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)
if (!length(files)) stop("No calibration result files found in ", cli$input)
raw <- do.call(rbind, lapply(files, read_tsv))
raw <- raw[raw$converged %in% TRUE & raw$he_converged %in% TRUE &
               is.finite(raw$true_h2) & is.finite(raw$he_h2), , drop = FALSE]
if (!nrow(raw)) stop("No successful calibration results")
if (!"null_alpha" %in% names(raw)) {
    warning("Legacy raw files lack null_alpha; using 0.05")
    raw$null_alpha <- 0.05
}
null_alpha_values <- unique(raw$null_alpha)
if (length(null_alpha_values) != 1L || !is.finite(null_alpha_values)) {
    stop("Calibration files must contain one finite null_alpha")
}
null_alpha <- null_alpha_values[[1L]]
raw_metrics <- unique(raw$raw_metric)
if (length(raw_metrics) != 1L) stop("Calibration files contain multiple raw metrics")
raw_metric <- raw_metrics[[1L]]
if (!raw_metric %in% names(raw)) stop("Raw metric column is absent: ", raw_metric)

setting_columns <- c("outer_folds", "outer_repeats", "inner_folds",
                     "alpha_grid", "lambda_rule", "max_features")
for (column in setting_columns) {
    if (length(unique(raw[[column]])) != 1L) {
        stop("Estimator setting varies across calibration results: ", column)
    }
}

criteria <- read_tsv(cli$criteria)
gate_versions <- unique(criteria$gate_version)
if (length(gate_versions) != 1L || !nzchar(gate_versions[[1L]])) {
    stop("Acceptance criteria must declare exactly one gate_version")
}
gate_version <- gate_versions[[1L]]
analysis_config <- read_tsv(cli$analysis_config)
analysis_settings <- stats::setNames(analysis_config$value, analysis_config$setting)
if (!"low_h2_max" %in% names(analysis_settings)) {
    stop("Analysis configuration must lock low_h2_max")
}
low_h2_max <- as_num(analysis_settings[["low_h2_max"]], "low_h2_max")
criterion <- function(metric) {
    value <- criteria$threshold[criteria$metric == metric]
    if (length(value) != 1L || !is.finite(value)) {
        stop("Calibration criterion is missing or invalid: ", metric)
    }
    as.numeric(value)
}
limits <- list(
    null_mean = criterion("null_mean_estimated_h2"),
    low_mean_bias = criterion("low_h2_absolute_mean_bias"),
    low_level_bias = criterion("max_absolute_low_h2_level_bias"),
    level_bias = criterion("max_absolute_h2_level_bias"),
    ordering = criterion("spearman_truth_estimate"),
    rmse_guardrail = criterion("rmse")
)
calibration_upper_bound <- max(raw$true_h2)

## Use a deterministic, design-balanced internal split. The evaluation split is
## never inspected while choosing the ensemble weight.
if ("replicate" %in% names(raw) && length(unique(raw$replicate)) >= 2L) {
    tune_index <- raw$replicate %% 2L == 0L
} else {
    design_key <- interaction(
        raw$n, raw$num_snps, raw$ld_rho, raw$architecture, raw$true_h2,
        drop = TRUE, lex.order = TRUE
    )
    tune_index <- ave(seq_len(nrow(raw)), design_key, FUN = seq_along) %% 2L == 0L
}
if (sum(tune_index) < 3L || sum(!tune_index) < 3L) {
    stop("Calibration results are too sparse for internal fit/tuning separation")
}
fit_data <- raw[!tune_index, , drop = FALSE]
tune_data <- raw[tune_index, , drop = FALSE]

internal_forward <- fit_forward_regression(fit_data)
fit_score <- predict_forward_regression(internal_forward, fit_data)$estimate
internal_debias <- fit_affine_level_debiasing(fit_score, fit_data$true_h2)
tune_score <- predict_forward_regression(internal_forward, tune_data)$estimate
tune_forward_h2 <- predict_affine_level_debiasing(internal_debias, tune_score)

weight_grid <- seq(0, 1, by = 0.025)
tuning <- do.call(rbind, lapply(weight_grid, function(weight) {
    estimate <- pmin(
        calibration_upper_bound,
        (1 - weight) * tune_forward_h2 + weight * tune_data$he_h2
    )
    error <- estimate - tune_data$true_h2
    level_bias <- aggregate(error, by = list(true_h2 = tune_data$true_h2), FUN = mean)
    low <- tune_data$true_h2 <= low_h2_max
    low_level_bias <- level_bias[level_bias$true_h2 <= low_h2_max, , drop = FALSE]
    null_estimate <- estimate[tune_data$true_h2 == 0]
    null_mean <- mean(null_estimate)
    null_mean_se <- stats::sd(null_estimate) / sqrt(length(null_estimate))
    data.frame(
        he_weight = weight,
        forward_weight = 1 - weight,
        null_mean_estimated_h2 = null_mean,
        null_mean_standard_error = null_mean_se,
        null_mean_upper_95 = null_mean + stats::qnorm(0.95) * null_mean_se,
        low_h2_absolute_mean_bias = abs(mean(error[low])),
        max_absolute_low_h2_level_bias = max(abs(low_level_bias$x)),
        max_absolute_h2_level_bias = max(abs(level_bias$x)),
        rmse = sqrt(mean(error^2)),
        spearman_truth_estimate = suppressWarnings(stats::cor(
            estimate, tune_data$true_h2, method = "spearman"
        )),
        stringsAsFactors = FALSE
    )
}))
tuning$rmse_guardrail_met <- tuning$rmse <= limits$rmse_guardrail
tuning$passes_internal_constraints <-
    tuning$null_mean_upper_95 <= limits$null_mean &
    tuning$low_h2_absolute_mean_bias <= limits$low_mean_bias &
    tuning$max_absolute_low_h2_level_bias <= limits$low_level_bias &
    tuning$max_absolute_h2_level_bias <= limits$level_bias &
    tuning$spearman_truth_estimate >= limits$ordering
if (any(tuning$passes_internal_constraints)) {
    eligible <- which(tuning$passes_internal_constraints)
    ## RMSE remains a useful optimization target among scientifically eligible
    ## weights, but crossing its guardrail cannot by itself reject a model.
    selected_row <- eligible[[which.min(tuning$rmse[eligible])]]
    selection_status <- "constraints_satisfied"
} else {
    ## Keep small smoke grids executable, but record that their internal tuning
    ## sample was insufficient to certify the locked full-analysis criteria.
    violation <- pmax(0, tuning$null_mean_upper_95 / limits$null_mean - 1) +
        pmax(0, tuning$low_h2_absolute_mean_bias / limits$low_mean_bias - 1) +
        pmax(0, tuning$max_absolute_low_h2_level_bias / limits$low_level_bias - 1) +
        pmax(0, tuning$max_absolute_h2_level_bias / limits$level_bias - 1) +
        pmax(0, limits$ordering / tuning$spearman_truth_estimate - 1)
    selected_row <- which.min(violation + tuning$rmse)
    selection_status <- "no_feasible_weight"
    warning("No hybrid weight met all internal constraints; final evaluation must reject this model")
}
tuning$selected <- seq_len(nrow(tuning)) == selected_row
selected_weight <- tuning$he_weight[[selected_row]]

## The interval residuals remain honestly out of sample with respect to the
## internal forward model and selected blend.
tune_estimate <- pmin(
    calibration_upper_bound,
    (1 - selected_weight) * tune_forward_h2 + selected_weight * tune_data$he_h2
)
tune_residual <- tune_data$true_h2 - tune_estimate
residual_quantiles <- stats::quantile(
    tune_residual, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE
)

## Refit the point-estimation components on all calibration simulations after
## choosing the weight. The untouched evaluation simulations remain independent.
final_forward <- fit_forward_regression(raw)
full_score <- predict_forward_regression(final_forward, raw)$estimate
final_debias <- fit_affine_level_debiasing(full_score, raw$true_h2)

raw$stratum_id <- forward_design_id(raw$n, raw$num_snps, raw$ld_rho)
split_data <- split(raw, raw$stratum_id)
strata <- lapply(names(split_data), function(id) {
    data <- split_data[[id]]
    null_values <- data[[raw_metric]][data$true_h2 == 0]
    if (length(null_values) < 3L) stop("Too few null simulations in stratum ", id)
    null_cutoff <- finite_sample_upper_threshold(null_values, alpha = null_alpha)
    scaler <- final_forward$scalers[[id]]
    list(
        id = id,
        n_center = scaler$n_center,
        p_center = scaler$p_center,
        ld_center = scaler$ld_center,
        ld_rho_design = scaler$ld_rho_design,
        raw_range = range(data[[raw_metric]], na.rm = TRUE),
        null_alpha = null_alpha,
        null_threshold_method = "split_conformal_order_statistic",
        null_raw_threshold_95 = null_cutoff$threshold,
        null_order_index = null_cutoff$order_index,
        null_simulations = null_cutoff$n,
        null_attainable_alpha = null_cutoff$attainable_alpha,
        n_simulations = nrow(data),
        architectures = sort(unique(data$architecture))
    )
})
names(strata) <- vapply(strata, `[[`, character(1L), "id")

model <- list(
    calibration_version = "forward_hybrid_v1",
    method = "nested_cross_fitted_elastic_net_forward_calibration_with_he_control_variate",
    estimand = "simulation_calibrated_local_snp_explained_methylation_variance",
    raw_metric = raw_metric,
    null_alpha = null_alpha,
    max_design_distance = unique(raw$max_design_distance)[[1L]],
    estimator_settings = as.list(raw[1L, setting_columns, drop = FALSE]),
    simulation_architectures = sort(unique(raw$architecture)),
    simulation_h2_values = sort(unique(raw$true_h2)),
    created_utc = format(Sys.time(), tz = "UTC"),
    internal_split = "odd_replicates_fit_even_replicates_tune",
    internal_selection_status = selection_status,
    acceptance_gate_version = gate_version,
    low_h2_max = low_h2_max,
    he_weight = selected_weight,
    forward_weight = 1 - selected_weight,
    calibration_upper_bound = calibration_upper_bound,
    forward_model = final_forward,
    affine_debiasing = final_debias,
    residual_quantiles = residual_quantiles,
    tuning_performance = tuning[selected_row, , drop = FALSE],
    strata = strata
)

dir.create(dirname(cli$output_model), recursive = TRUE, showWarnings = FALSE)
tmp_model <- tempfile(pattern = ".calibration-", tmpdir = dirname(cli$output_model), fileext = ".rds")
saveRDS(model, tmp_model, version = 3)
if (!file.rename(tmp_model, cli$output_model)) stop("Could not write calibration model")

manifest <- do.call(rbind, lapply(strata, function(x) data.frame(
    calibration_stratum = x$id,
    n_center = x$n_center,
    p_center = x$p_center,
    ld_metric_center = x$ld_center,
    ld_rho_design = x$ld_rho_design,
    raw_metric_min = x$raw_range[[1L]],
    raw_metric_max = x$raw_range[[2L]],
    null_raw_threshold_95 = x$null_raw_threshold_95,
    null_alpha = x$null_alpha,
    null_threshold_method = x$null_threshold_method,
    null_order_index = x$null_order_index,
    null_simulations = x$null_simulations,
    null_attainable_alpha = x$null_attainable_alpha,
    residual_q025 = residual_quantiles[[1L]],
    residual_q975 = residual_quantiles[[2L]],
    he_weight = selected_weight,
    forward_weight = 1 - selected_weight,
    calibration_upper_bound = calibration_upper_bound,
    internal_selection_status = selection_status,
    acceptance_gate_version = gate_version,
    low_h2_max = low_h2_max,
    rmse_guardrail_met = tuning$rmse_guardrail_met[[selected_row]],
    n_simulations = x$n_simulations,
    architectures = paste(x$architectures, collapse = ","),
    stringsAsFactors = FALSE
)))
write_tsv(manifest, cli$output_manifest)
write_tsv(tuning, cli$output_tuning)
capture_session_info(cli$session_info)
cat("Wrote hybrid forward calibration model:", normalizePath(cli$output_model), "\n")
cat("Selected HE weight:", selected_weight, "(", selection_status, ")\n")
