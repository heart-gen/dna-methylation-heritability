#!/usr/bin/env Rscript

## Calibration-only model development. This script never reads raw/evaluation.
## It compares prespecified nonlinear forms on the odd-replicate fit / even-
## replicate tuning split and writes diagnostics; independent evaluation is a
## separate, one-time gate.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    input = "",
    criteria = file.path(dirname(script_path), "..", "config",
                         "acceptance-criteria.tsv"),
    analysis_config = file.path(dirname(script_path), "..", "config",
                                "analysis.tsv"),
    output = ""
))
if (!nzchar(cli$input) || !nzchar(cli$output)) {
    stop("--input and --output are required")
}
raw <- read_tsv(cli$input)
raw <- raw[raw$converged %in% TRUE & raw$he_converged %in% TRUE, , drop = FALSE]
fit_data <- raw[raw$replicate %% 2L == 1L, , drop = FALSE]
tune_data <- raw[raw$replicate %% 2L == 0L, , drop = FALSE]
criteria <- read_tsv(cli$criteria)
threshold <- stats::setNames(as.numeric(criteria$threshold), criteria$metric)
analysis_config <- read_tsv(cli$analysis_config)
analysis_settings <- stats::setNames(analysis_config$value, analysis_config$setting)
if (!"low_h2_max" %in% names(analysis_settings)) {
    stop("Analysis configuration must lock low_h2_max")
}
low_h2_max <- as_num(analysis_settings[["low_h2_max"]], "low_h2_max")
upper_bound <- max(raw$true_h2)

features_for_model <- function(data, replacements = NULL) {
    features <- make_forward_features(data)
    if (is.null(replacements)) {
        replacements <- vapply(features, function(x) {
            value <- stats::median(x[is.finite(x)])
            if (is.finite(value)) value else 0
        }, numeric(1L))
    }
    for (column in names(features)) {
        features[[column]][!is.finite(features[[column]])] <- replacements[[column]]
    }
    features$true_h2 <- data$true_h2
    features$stratum_id <- forward_design_id(data$n, data$num_snps, data$ld_rho)
    attr(features, "replacements") <- replacements
    features
}
fit_features <- features_for_model(fit_data)
tune_features <- features_for_model(
    tune_data, replacements = attr(fit_features, "replacements")
)

fit_per_stratum <- function(degree) {
    split_fit <- split(seq_len(nrow(fit_features)), fit_features$stratum_id)
    polynomial_terms <- paste0(
        "I(sqrt_rho^", seq_len(degree), ")", collapse = " + "
    )
    formula <- stats::as.formula(paste(
        "true_h2 ~", polynomial_terms,
        "+ r2 + covariance + sqrt_score_variance + log_nonzero_snps"
    ))
    models <- lapply(split_fit, function(index) {
        stats::lm(formula, data = fit_features[index, , drop = FALSE])
    })
    pred_fit <- numeric(nrow(fit_features))
    pred_tune <- numeric(nrow(tune_features))
    for (id in names(models)) {
        fit_index <- which(fit_features$stratum_id == id)
        tune_index <- which(tune_features$stratum_id == id)
        pred_fit[fit_index] <- stats::predict(
            models[[id]], newdata = fit_features[fit_index, , drop = FALSE]
        )
        pred_tune[tune_index] <- stats::predict(
            models[[id]], newdata = tune_features[tune_index, , drop = FALSE]
        )
    }
    debias <- fit_affine_level_debiasing(pred_fit, fit_data$true_h2)
    list(
        fit = predict_affine_level_debiasing(debias, pred_fit),
        tune = predict_affine_level_debiasing(debias, pred_tune)
    )
}

fit_pooled_spline <- function(df) {
    fit_features$n_factor <- factor(fit_data$n)
    fit_features$p_factor <- factor(fit_data$num_snps)
    fit_features$ld_factor <- factor(fit_data$ld_rho)
    tune_features$n_factor <- factor(tune_data$n, levels = levels(fit_features$n_factor))
    tune_features$p_factor <- factor(tune_data$num_snps, levels = levels(fit_features$p_factor))
    tune_features$ld_factor <- factor(tune_data$ld_rho, levels = levels(fit_features$ld_factor))
    formula <- true_h2 ~
        (splines::ns(sqrt_rho, df = df) + splines::ns(r2, df = df) +
         splines::ns(covariance, df = df) +
         splines::ns(sqrt_score_variance, df = df) +
         splines::ns(log_nonzero_snps, df = df)) *
        (n_factor + p_factor + ld_factor)
    model <- stats::lm(formula, data = fit_features)
    pred_fit <- as.numeric(stats::predict(model, newdata = fit_features))
    pred_tune <- as.numeric(stats::predict(model, newdata = tune_features))
    debias <- fit_affine_level_debiasing(pred_fit, fit_data$true_h2)
    list(
        fit = predict_affine_level_debiasing(debias, pred_fit),
        tune = predict_affine_level_debiasing(debias, pred_tune)
    )
}

score_candidate <- function(name, predictions) {
    do.call(rbind, lapply(seq(0, 1, by = 0.025), function(weight) {
        estimate <- pmin(
            upper_bound,
            (1 - weight) * predictions$tune + weight * tune_data$he_h2
        )
        error <- estimate - tune_data$true_h2
        level <- aggregate(error, by = list(true_h2 = tune_data$true_h2), FUN = mean)
        low <- tune_data$true_h2 <= low_h2_max
        low_level <- level[level$true_h2 <= low_h2_max, , drop = FALSE]
        null <- estimate[tune_data$true_h2 == 0]
        null_se <- stats::sd(null) / sqrt(length(null))
        data.frame(
            candidate = name,
            he_weight = weight,
            null_mean_estimated_h2 = mean(null),
            null_mean_upper_95 = mean(null) + stats::qnorm(0.95) * null_se,
            low_h2_absolute_mean_bias = abs(mean(error[low])),
            max_absolute_low_h2_level_bias = max(abs(low_level$x)),
            max_absolute_h2_level_bias = max(abs(level$x)),
            rmse = sqrt(mean(error^2)),
            absolute_mean_bias = abs(mean(error)),
            spearman_truth_estimate = suppressWarnings(stats::cor(
                estimate, tune_data$true_h2, method = "spearman"
            )),
            stringsAsFactors = FALSE
        )
    }))
}

results <- list()
index <- 1L
for (df in 1:4) {
    results[[index]] <- score_candidate(
        paste0("per_stratum_polynomial_degree", df), fit_per_stratum(df)
    )
    index <- index + 1L
}
for (df in 2:5) {
    results[[index]] <- score_candidate(
        paste0("pooled_interaction_spline_df", df), fit_pooled_spline(df)
    )
    index <- index + 1L
}
results <- do.call(rbind, results)
results$passes_internal_constraints <-
    results$null_mean_upper_95 <= threshold[["null_mean_estimated_h2"]] &
    results$low_h2_absolute_mean_bias <=
        threshold[["low_h2_absolute_mean_bias"]] &
    results$max_absolute_low_h2_level_bias <=
        threshold[["max_absolute_low_h2_level_bias"]] &
    results$max_absolute_h2_level_bias <=
        threshold[["max_absolute_h2_level_bias"]] &
    results$spearman_truth_estimate >= threshold[["spearman_truth_estimate"]]
results$rmse_guardrail_met <- results$rmse <= threshold[["rmse"]]
results$selected <- FALSE
eligible <- which(results$passes_internal_constraints)
if (length(eligible)) {
    selected <- eligible[[which.min(results$rmse[eligible])]]
    results$selected[[selected]] <- TRUE
}
write_tsv(results, cli$output)
cat("Internal candidates passing:", length(eligible), "of", nrow(results), "\n")
if (length(eligible)) print(results[results$selected, , drop = FALSE])
