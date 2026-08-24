#!/usr/bin/env Rscript

## Fit prespecified boundary-aware / two-part calibrators on a calibration-only
## training set. Selection uses the odd-fit / even-tune split only. The
## motivating failed validation grid is never read.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    input = file.path(dirname(script_path), "..", "_m", "raw", "calibration"),
    candidates = file.path(dirname(script_path), "..", "config",
                          "boundary-aware-candidates.tsv"),
    criteria = file.path(dirname(script_path), "..", "config",
                         "acceptance-criteria.tsv"),
    analysis_config = file.path(dirname(script_path), "..", "config",
                                "analysis.tsv"),
    output_model = file.path(dirname(script_path), "..", "_m", "calibration",
                             "elastic-net-calibration.rds"),
    output_manifest = file.path(dirname(script_path), "..", "_m", "calibration",
                                "calibration-manifest.tsv"),
    output_tuning = file.path(dirname(script_path), "..", "_m", "calibration",
                              "candidate-selection.tsv"),
    output_selected = file.path(dirname(script_path), "..", "_m", "calibration",
                                "selected-candidate.tsv"),
    session_info = file.path(dirname(script_path), "..", "_m", "calibration",
                             "session-info.txt"),
    fail_closed = "TRUE"
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
max_simulated_h2 <- max(raw$true_h2)
candidates <- read_tsv(cli$candidates)
required_candidate_cols <- c(
    "candidate_id", "family", "parameter", "upper_bound_rule"
)
missing_candidate <- setdiff(required_candidate_cols, names(candidates))
if (length(missing_candidate)) {
    stop("Candidate menu missing columns: ", paste(missing_candidate, collapse = ", "))
}

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
weight_grid <- seq(0, 1, by = 0.025)

tuning_rows <- list()
for (i in seq_len(nrow(candidates))) {
    cand <- candidates[i, , drop = FALSE]
    upper_bound <- resolve_calibration_upper_bound(
        cand$upper_bound_rule[[1L]], max_simulated_h2
    )
    candidate_fit <- tryCatch(
        fit_candidate_forward_score(
            cand$candidate_id[[1L]],
            cand$family[[1L]],
            cand$parameter[[1L]],
            fit_data
        ),
        error = function(e) e
    )
    if (inherits(candidate_fit, "error")) {
        tuning_rows[[length(tuning_rows) + 1L]] <- data.frame(
            candidate_id = cand$candidate_id[[1L]],
            family = cand$family[[1L]],
            parameter = as.character(cand$parameter[[1L]]),
            upper_bound_rule = cand$upper_bound_rule[[1L]],
            he_weight = NA_real_,
            null_mean_estimated_h2 = NA_real_,
            null_mean_standard_error = NA_real_,
            null_mean_upper_95 = NA_real_,
            low_h2_absolute_mean_bias = NA_real_,
            max_absolute_low_h2_level_bias = NA_real_,
            max_absolute_h2_level_bias = NA_real_,
            rmse = NA_real_,
            absolute_mean_bias = NA_real_,
            spearman_truth_estimate = NA_real_,
            rmse_guardrail_met = FALSE,
            passes_internal_constraints = FALSE,
            selected = FALSE,
            fit_error = conditionMessage(candidate_fit),
            stringsAsFactors = FALSE
        )
        next
    }
    fit_score <- predict_candidate_forward_score(candidate_fit, fit_data)$estimate
    debias <- fit_affine_level_debiasing(fit_score, fit_data$true_h2)
    tune_score <- predict_candidate_forward_score(candidate_fit, tune_data)$estimate
    tune_forward <- predict_affine_level_debiasing(debias, tune_score)
    scored <- do.call(rbind, lapply(weight_grid, function(weight) {
        row <- score_hybrid_on_tune(
            tune_forward, tune_data$he_h2, tune_data$true_h2, upper_bound, weight,
            low_h2_max = low_h2_max
        )
        row$candidate_id <- cand$candidate_id[[1L]]
        row$family <- cand$family[[1L]]
        row$parameter <- as.character(cand$parameter[[1L]])
        row$upper_bound_rule <- cand$upper_bound_rule[[1L]]
        row$fit_error <- NA_character_
        row
    }))
    scored$rmse_guardrail_met <- scored$rmse <= limits$rmse_guardrail
    scored$passes_internal_constraints <-
        is.finite(scored$null_mean_upper_95) &
        scored$null_mean_upper_95 <= limits$null_mean &
        scored$low_h2_absolute_mean_bias <= limits$low_mean_bias &
        scored$max_absolute_low_h2_level_bias <= limits$low_level_bias &
        scored$max_absolute_h2_level_bias <= limits$level_bias &
        scored$spearman_truth_estimate >= limits$ordering
    scored$selected <- FALSE
    tuning_rows[[length(tuning_rows) + 1L]] <- scored[
        ,
        c(
            "candidate_id", "family", "parameter", "upper_bound_rule", "he_weight",
            "null_mean_estimated_h2", "null_mean_standard_error",
            "null_mean_upper_95", "low_h2_absolute_mean_bias",
            "max_absolute_low_h2_level_bias", "max_absolute_h2_level_bias", "rmse",
            "absolute_mean_bias", "spearman_truth_estimate",
            "rmse_guardrail_met", "passes_internal_constraints", "selected", "fit_error"
        ),
        drop = FALSE
    ]
}
tuning <- do.call(rbind, tuning_rows)
eligible <- which(tuning$passes_internal_constraints %in% TRUE)
if (!length(eligible)) {
    write_tsv(tuning, cli$output_tuning)
    capture_session_info(cli$session_info)
    msg <- paste0(
        "No boundary-aware/two-part candidate met internal constraints on the ",
        "training tune split; refusing to freeze a calibrator"
    )
    if (identical(toupper(cli$fail_closed), "TRUE")) {
        stop(msg)
    }
    warning(msg)
    quit(save = "no", status = 2)
}
selected_row <- eligible[[which.min(tuning$rmse[eligible])]]
tuning$selected[[selected_row]] <- TRUE
selected <- tuning[selected_row, , drop = FALSE]
selected$acceptance_gate_version <- gate_version
selected$low_h2_max <- low_h2_max
selected_weight <- selected$he_weight[[1L]]
selected_upper <- resolve_calibration_upper_bound(
    selected$upper_bound_rule[[1L]], max_simulated_h2
)

## Refit selected family on all calibration simulations after weight choice.
final_candidate <- fit_candidate_forward_score(
    selected$candidate_id[[1L]],
    selected$family[[1L]],
    selected$parameter[[1L]],
    raw
)
full_score <- predict_candidate_forward_score(final_candidate, raw)$estimate
final_debias <- fit_affine_level_debiasing(full_score, raw$true_h2)

## Honest out-of-sample residuals from the tune split of the selected form.
selected_fit <- fit_candidate_forward_score(
    selected$candidate_id[[1L]],
    selected$family[[1L]],
    selected$parameter[[1L]],
    fit_data
)
fit_score <- predict_candidate_forward_score(selected_fit, fit_data)$estimate
tune_debias <- fit_affine_level_debiasing(fit_score, fit_data$true_h2)
tune_score <- predict_candidate_forward_score(selected_fit, tune_data)$estimate
tune_forward <- predict_affine_level_debiasing(tune_debias, tune_score)
tune_estimate <- pmin(
    selected_upper,
    (1 - selected_weight) * tune_forward + selected_weight * tune_data$he_h2
)
tune_residual <- tune_data$true_h2 - tune_estimate
residual_quantiles <- stats::quantile(
    tune_residual, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE
)

raw$stratum_id <- forward_design_id(raw$n, raw$num_snps, raw$ld_rho)
## Prefer strata centers from the selected family's primary forward scalers.
primary_scalers <- if (!is.null(final_candidate$forward_model)) {
    final_candidate$forward_model$scalers
} else if (!is.null(final_candidate$interior_model)) {
    final_candidate$interior_model$scalers
} else if (!is.null(final_candidate$low_model)) {
    final_candidate$low_model$scalers
} else {
    fit_forward_scalers(raw)
}
split_data <- split(raw, raw$stratum_id)
strata <- lapply(names(split_data), function(id) {
    data <- split_data[[id]]
    null_values <- data[[raw_metric]][data$true_h2 == 0]
    if (length(null_values) < 3L) stop("Too few null simulations in stratum ", id)
    null_cutoff <- finite_sample_upper_threshold(null_values, alpha = null_alpha)
    scaler <- primary_scalers[[id]]
    if (is.null(scaler)) {
        scaler <- list(
            n_center = stats::median(data$n),
            p_center = stats::median(data$num_snps),
            ld_center = stats::median(data$ld_metric, na.rm = TRUE),
            ld_rho_design = unique(data$ld_rho)[[1L]]
        )
    }
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
    calibration_version = "candidate_hybrid_v1",
    candidate_id = selected$candidate_id[[1L]],
    candidate_family = selected$family[[1L]],
    candidate_parameter = selected$parameter[[1L]],
    method = "nested_cross_fitted_elastic_net_boundary_aware_candidate_calibration",
    estimand = "simulation_calibrated_local_snp_explained_methylation_variance",
    raw_metric = raw_metric,
    null_alpha = null_alpha,
    max_design_distance = unique(raw$max_design_distance)[[1L]],
    estimator_settings = as.list(raw[1L, setting_columns, drop = FALSE]),
    simulation_architectures = sort(unique(raw$architecture)),
    simulation_h2_values = sort(unique(raw$true_h2)),
    created_utc = format(Sys.time(), tz = "UTC"),
    internal_split = "odd_replicates_fit_even_replicates_tune",
    internal_selection_status = "constraints_satisfied",
    acceptance_gate_version = gate_version,
    low_h2_max = low_h2_max,
    he_weight = selected_weight,
    forward_weight = 1 - selected_weight,
    calibration_upper_bound = selected_upper,
    candidate_model = final_candidate,
    affine_debiasing = final_debias,
    residual_quantiles = residual_quantiles,
    tuning_performance = selected,
    strata = strata
)

dir.create(dirname(cli$output_model), recursive = TRUE, showWarnings = FALSE)
tmp_model <- tempfile(
    pattern = ".calibration-", tmpdir = dirname(cli$output_model), fileext = ".rds"
)
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
    calibration_upper_bound = selected_upper,
    candidate_id = selected$candidate_id[[1L]],
    candidate_family = selected$family[[1L]],
    internal_selection_status = "constraints_satisfied",
    acceptance_gate_version = gate_version,
    low_h2_max = low_h2_max,
    rmse_guardrail_met = selected$rmse_guardrail_met[[1L]],
    n_simulations = x$n_simulations,
    architectures = paste(x$architectures, collapse = ","),
    stringsAsFactors = FALSE
)))
write_tsv(manifest, cli$output_manifest)
write_tsv(tuning, cli$output_tuning)
write_tsv(selected, cli$output_selected)
capture_session_info(cli$session_info)
cat(
    "Froze candidate", selected$candidate_id[[1L]],
    "HE weight", selected_weight,
    "upper bound", selected_upper, "\n"
)
cat("Wrote model:", normalizePath(cli$output_model), "\n")
