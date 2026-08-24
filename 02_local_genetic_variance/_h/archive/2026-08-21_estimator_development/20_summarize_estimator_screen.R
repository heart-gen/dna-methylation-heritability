#!/usr/bin/env Rscript

## Estimator-screening summary. Compares settings arms on identical simulated
## draws. This never writes a calibration model and never touches production.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    input = "",
    config = file.path(dirname(script_path), "..", "config",
                       "estimator-screen-20260818.tsv"),
    manifest = "",
    output_dir = "",
    fail_on_rejection = "FALSE"
))
if (!nzchar(cli$input) || !nzchar(cli$output_dir) || !nzchar(cli$manifest)) {
    stop("--input, --manifest, and --output-dir are required")
}
config_table <- read_tsv(cli$config)
settings <- stats::setNames(as.list(config_table$value), config_table$field)
improvement_rule <- 0.15
reference_arm <- "A0_baseline"

files <- list.files(cli$input, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)
if (!length(files)) stop("No screening result files found in ", cli$input)
raw <- do.call(rbind, lapply(files, read_tsv))
manifest <- read_tsv(cli$manifest)
raw <- merge(
    raw, manifest[, c("scenario_id", "arm_id", "cell_id")],
    by = "scenario_id", all.x = TRUE
)

dir.create(cli$output_dir, recursive = TRUE, showWarnings = FALSE)

## Reconciliation before any scientific readout.
reconciliation <- data.frame(
    expected = nrow(manifest),
    observed = nrow(raw),
    converged = sum(raw$converged %in% TRUE),
    he_converged = sum(raw$he_converged %in% TRUE),
    missing = nrow(manifest) - nrow(raw),
    stringsAsFactors = FALSE
)
write_tsv(reconciliation, file.path(cli$output_dir, "screen-reconciliation.tsv"))

## Negative control: Haseman-Elston is settings-invariant, so matched cells must
## agree across arms. Disagreement means the pairing is broken.
he_check <- aggregate(he_h2 ~ cell_id, raw, function(x) {
    if (all(is.na(x))) 0 else diff(range(x, na.rm = TRUE))
})
names(he_check)[[2L]] <- "he_range_across_arms"
pairing_max_range <- max(he_check$he_range_across_arms, na.rm = TRUE)
pairing_ok <- is.finite(pairing_max_range) && pairing_max_range < 1e-8
write_tsv(he_check, file.path(cli$output_dir, "pairing-negative-control.tsv"))
cat("Pairing negative control: max HE range across arms =",
    signif(pairing_max_range, 3), if (pairing_ok) "(OK)" else "(BROKEN)", "\n")

analysis <- raw[raw$converged %in% TRUE & raw$he_converged %in% TRUE &
                    is.finite(raw$true_h2) & is.finite(raw$he_h2), , drop = FALSE]

## Secondary endpoint: resolvability of the raw statistic. A calibration map
## cannot separate h2 levels that the raw metric does not separate.
discriminability <- do.call(rbind, lapply(
    split(analysis, list(analysis$arm_id, analysis$n, analysis$num_snps,
                         analysis$architecture), drop = TRUE),
    function(d) {
        level_mean <- tapply(d$rho2_oof, d$true_h2, mean, na.rm = TRUE)
        level_sd <- tapply(d$rho2_oof, d$true_h2, stats::sd, na.rm = TRUE)
        within <- mean(level_sd, na.rm = TRUE)
        between <- stats::sd(level_mean, na.rm = TRUE)
        data.frame(
            arm_id = d$arm_id[[1L]], n = d$n[[1L]], num_snps = d$num_snps[[1L]],
            architecture = d$architecture[[1L]],
            within_level_sd = within, between_level_sd = between,
            discriminability = safe_ratio(between, within),
            simulations = nrow(d), stringsAsFactors = FALSE
        )
    }
))
rownames(discriminability) <- NULL
write_tsv(discriminability, file.path(cli$output_dir, "screen-discriminability.tsv"))

## Primary endpoint: oracle calibrated RMSE. Within each arm, fit the forward
## map on odd replicates and score even replicates, mirroring the production
## internal split. Forward-only isolates the elastic-net component that the
## arms actually change; the fixed-weight hybrid is reported alongside because
## production reports a hybrid.
hybrid_weight <- 0.175
score_arm <- function(d) {
    fit_data <- d[d$replicate %% 2L == 1L, , drop = FALSE]
    tune_data <- d[d$replicate %% 2L == 0L, , drop = FALSE]
    if (nrow(fit_data) < 30L || nrow(tune_data) < 30L) return(NULL)
    model <- tryCatch(fit_forward_regression(fit_data), error = function(e) e)
    if (inherits(model, "error")) return(NULL)
    fit_score <- predict_forward_regression(model, fit_data)$estimate
    debias <- tryCatch(fit_affine_level_debiasing(fit_score, fit_data$true_h2),
                       error = function(e) e)
    if (inherits(debias, "error")) return(NULL)
    tune_score <- predict_forward_regression(model, tune_data)$estimate
    forward <- predict_affine_level_debiasing(debias, tune_score)
    upper <- max(d$true_h2)
    forward_estimate <- pmin(upper, forward)
    hybrid_estimate <- pmin(
        upper, (1 - hybrid_weight) * forward + hybrid_weight * tune_data$he_h2
    )
    data.frame(
        arm_id = d$arm_id[[1L]],
        scenario_id = tune_data$scenario_id,
        cell_id = tune_data$cell_id,
        n = tune_data$n, num_snps = tune_data$num_snps,
        architecture = tune_data$architecture, true_h2 = tune_data$true_h2,
        rho2_oof = tune_data$rho2_oof, he_h2 = tune_data$he_h2,
        forward_estimate = forward_estimate,
        hybrid_estimate = hybrid_estimate,
        forward_error = forward_estimate - tune_data$true_h2,
        hybrid_error = hybrid_estimate - tune_data$true_h2,
        stringsAsFactors = FALSE
    )
}
scored <- do.call(rbind, lapply(split(analysis, analysis$arm_id), score_arm))
if (is.null(scored) || !nrow(scored)) stop("No arm could be scored")
rownames(scored) <- NULL
write_tsv(scored, file.path(cli$output_dir, "screen-scored-estimates.tsv"))

summarize_arm <- function(d) {
    null_rows <- d[d$true_h2 == 0, , drop = FALSE]
    level_bias <- tapply(d$forward_error, d$true_h2, mean, na.rm = TRUE)
    data.frame(
        arm_id = d$arm_id[[1L]],
        simulations = nrow(d),
        oracle_rmse_forward = sqrt(mean(d$forward_error^2, na.rm = TRUE)),
        oracle_rmse_hybrid = sqrt(mean(d$hybrid_error^2, na.rm = TRUE)),
        absolute_mean_bias = abs(mean(d$forward_error, na.rm = TRUE)),
        max_absolute_h2_level_bias = max(abs(level_bias), na.rm = TRUE),
        null_mean_estimated_h2 = if (nrow(null_rows)) {
            mean(null_rows$forward_estimate, na.rm = TRUE)
        } else {
            NA_real_
        },
        spearman_truth_estimate = suppressWarnings(stats::cor(
            d$forward_estimate, d$true_h2, method = "spearman",
            use = "complete.obs"
        )),
        stringsAsFactors = FALSE
    )
}
overall <- do.call(rbind, lapply(split(scored, scored$arm_id), summarize_arm))
rownames(overall) <- NULL

stress <- scored[scored$architecture == "polygenic" &
                     scored$num_snps == max(scored$num_snps) &
                     scored$true_h2 >= 0.6, , drop = FALSE]
stress_summary <- do.call(rbind, lapply(split(stress, stress$arm_id), function(d) {
    data.frame(
        arm_id = d$arm_id[[1L]], stress_simulations = nrow(d),
        stress_rmse = sqrt(mean(d$forward_error^2, na.rm = TRUE)),
        stress_bias = mean(d$forward_error, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
}))
rownames(stress_summary) <- NULL
overall <- merge(overall, stress_summary, by = "arm_id", all.x = TRUE)

if (!reference_arm %in% overall$arm_id) {
    stop("Reference arm is absent from the screen: ", reference_arm)
}
reference <- overall[overall$arm_id == reference_arm, , drop = FALSE]
overall$stress_rmse_relative_improvement <-
    (reference$stress_rmse - overall$stress_rmse) / reference$stress_rmse
overall$overall_rmse_relative_improvement <-
    (reference$oracle_rmse_forward - overall$oracle_rmse_forward) /
        reference$oracle_rmse_forward
overall$null_mean_worsened <-
    overall$null_mean_estimated_h2 > reference$null_mean_estimated_h2 + 0.01
overall$promising <- overall$arm_id != reference_arm &
    overall$stress_rmse_relative_improvement >= improvement_rule &
    !overall$null_mean_worsened
overall <- overall[order(overall$oracle_rmse_forward), ]
write_tsv(overall, file.path(cli$output_dir, "screen-arm-performance.tsv"))

promising <- overall$arm_id[overall$promising %in% TRUE]
decision <- data.frame(
    reference_arm = reference_arm,
    pairing_negative_control_passed = pairing_ok,
    reconciliation_complete = reconciliation$missing == 0L,
    improvement_rule = improvement_rule,
    promising_arms = if (length(promising)) {
        paste(promising, collapse = ",")
    } else {
        NA_character_
    },
    decision = if (!pairing_ok) {
        "INVALID_PAIRING_BROKEN"
    } else if (length(promising)) {
        "PROMISING_EXPAND_TO_RECALIBRATION_GRID"
    } else {
        "NO_TUNING_GAIN_INFORMATION_LIMIT"
    },
    note = paste(
        "Screening only. A promising arm still requires a full recalibration",
        "grid and a fresh validation split before any production change",
        "(AGENTS.md 7.2)."
    ),
    stringsAsFactors = FALSE
)
write_tsv(decision, file.path(cli$output_dir, "screen-decision.tsv"))
capture_session_info(file.path(cli$output_dir, "session-info.txt"))

print(overall[, c("arm_id", "oracle_rmse_forward", "stress_rmse",
                  "stress_rmse_relative_improvement",
                  "null_mean_estimated_h2", "promising")],
      row.names = FALSE, digits = 4)
cat("\nDecision:", decision$decision, "\n")
if (identical(toupper(cli$fail_on_rejection), "TRUE") &&
    !identical(decision$decision, "PROMISING_EXPAND_TO_RECALIBRATION_GRID")) {
    quit(save = "no", status = 2)
}
