#!/usr/bin/env Rscript
## Absolute Module 02 acceptance evaluation for raw BSLMM PVE.
## Does not modify observed-data production.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))
source(file.path(dirname(script_path), "bslmm_pilot_functions.R"))

cli <- parse_cli(list(
    input_dir = "",
    manifest = "",
    config = "",
    criteria = "",
    output_dir = "",
    fail_on_rejection = "TRUE"
))
if (!nzchar(cli$input_dir) || !nzchar(cli$manifest) ||
    !nzchar(cli$config) || !nzchar(cli$output_dir)) {
    stop("--input_dir, --manifest, --config, and --output_dir are required")
}
settings <- read_pilot_settings(cli$config)
manifest <- read_tsv(cli$manifest)
criteria_path <- if (nzchar(cli$criteria)) {
    cli$criteria
} else {
    file.path(dirname(script_path), "..", settings$acceptance_criteria_relpath)
}
## Allow criteria path relative to analysis dir
if (!file.exists(criteria_path)) {
    alt <- file.path(dirname(script_path), "..", basename(criteria_path))
    if (file.exists(alt)) criteria_path <- alt
}
if (!file.exists(criteria_path)) {
    ## config says config/acceptance-criteria.tsv relative to module
    criteria_path <- file.path(dirname(script_path), "..", "config",
                               "acceptance-criteria.tsv")
}
criteria <- read_tsv(criteria_path)
gate_versions <- unique(criteria$gate_version)
if (length(gate_versions) != 1L || !nzchar(gate_versions[[1L]])) {
    stop("Acceptance criteria must declare exactly one gate_version")
}
gate_version <- gate_versions[[1L]]

files <- list.files(cli$input_dir, pattern = "^scenario-[0-9]+\\.tsv$",
                    full.names = TRUE)
if (!length(files)) stop("No scenario files in ", cli$input_dir)
tables <- lapply(files, read_tsv)
cols <- Reduce(union, lapply(tables, names))
tables <- lapply(tables, function(d) {
    for (m in setdiff(cols, names(d))) d[[m]] <- NA
    d[, cols, drop = FALSE]
})
evaluation <- do.call(rbind, tables)
evaluation <- evaluation[order(evaluation$scenario_id), , drop = FALSE]
dir.create(cli$output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(evaluation, file.path(cli$output_dir, "bslmm-validation-estimates.tsv"))

missing_ids <- setdiff(manifest$scenario_id, evaluation$scenario_id)
reconciliation <- data.frame(
    expected_scenarios = nrow(manifest),
    completed_scenarios = nrow(evaluation),
    missing_scenarios = length(missing_ids),
    computational_failures = sum(evaluation$failed %in% TRUE),
    complete = length(missing_ids) == 0L,
    stringsAsFactors = FALSE
)
write_tsv(reconciliation, file.path(cli$output_dir, "task-reconciliation.tsv"))
if (length(missing_ids)) {
    write_tsv(data.frame(scenario_id = missing_ids),
              file.path(cli$output_dir, "missing-scenarios.tsv"))
}

## Map BSLMM columns onto the Module 02 acceptance metric names.
evaluation$h2_en_calibrated <- evaluation$bslmm_pve
evaluation$h2_calibration_lower <- evaluation$bslmm_pve_q025
evaluation$h2_calibration_upper <- evaluation$bslmm_pve_q975

mean_or_na <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
available <- is.finite(evaluation$bslmm_pve) & !(evaluation$failed %in% TRUE)
err <- evaluation$error
overall <- data.frame(
    simulations = nrow(evaluation),
    estimable_simulations = sum(available),
    estimable_rate = mean(available),
    mean_true_h2 = mean_or_na(evaluation$true_h2),
    mean_estimated_h2 = mean_or_na(evaluation$bslmm_pve),
    bias = mean_or_na(err),
    absolute_mean_bias = abs(mean_or_na(err)),
    mean_absolute_error = mean_or_na(abs(err)),
    rmse = sqrt(mean_or_na(err^2)),
    calibration_interval_coverage = mean_or_na(evaluation$covered),
    positive_signal_rate = mean_or_na(evaluation$positive_signal),
    within_domain_rate = mean(evaluation$within_domain %in% TRUE),
    computational_failure_rate = mean(evaluation$failed %in% TRUE),
    out_of_bounds_rate = mean_or_na(
        evaluation$bslmm_pve < 0 | evaluation$bslmm_pve > 1
    ),
    stringsAsFactors = FALSE
)
overall$spearman_truth_estimate <- suppressWarnings(stats::cor(
    evaluation$true_h2, evaluation$bslmm_pve,
    method = "spearman", use = "complete.obs"
))
null <- evaluation[evaluation$true_h2 == 0, , drop = FALSE]
overall$null_type1_error <- if (nrow(null)) {
    mean(null$positive_signal %in% TRUE)
} else {
    NA_real_
}
overall$null_mean_estimated_h2 <- if (nrow(null)) {
    mean_or_na(null$bslmm_pve)
} else {
    NA_real_
}
h2_level_bias <- aggregate(error ~ true_h2, evaluation, function(x) mean(x, na.rm = TRUE))
low_pve_max <- as_num(settings$low_pve_max %||% "0.1", "low_pve_max")
low <- evaluation$true_h2 <= low_pve_max
low_h2_level_bias <- h2_level_bias[
    h2_level_bias$true_h2 <= low_pve_max, , drop = FALSE
]
if (!any(low) || !nrow(low_h2_level_bias)) {
    stop("No BSLMM simulations fall in the configured low-h2 range")
}
overall$low_h2_max <- low_pve_max
overall$low_h2_absolute_mean_bias <- abs(mean_or_na(evaluation$error[low]))
overall$max_absolute_low_h2_level_bias <- max(
    abs(low_h2_level_bias$error), na.rm = TRUE
)
overall$max_absolute_h2_level_bias <- max(abs(h2_level_bias$error), na.rm = TRUE)
write_tsv(overall, file.path(cli$output_dir, "bslmm-validation-performance-overall.tsv"))
write_tsv(h2_level_bias, file.path(cli$output_dir, "bslmm-validation-bias-by-h2.tsv"))

group_columns <- c("n", "num_snps", "ld_rho", "architecture", "true_h2")
group_key <- interaction(evaluation[group_columns], drop = TRUE, lex.order = TRUE)
by_design <- do.call(rbind, lapply(split(evaluation, group_key), function(data) {
    avail <- is.finite(data$bslmm_pve) & !(data$failed %in% TRUE)
    e <- data$error
    cbind(
        data[1L, group_columns, drop = FALSE],
        data.frame(
            simulations = nrow(data),
            estimable_simulations = sum(avail),
            bias = mean_or_na(e),
            absolute_mean_bias = abs(mean_or_na(e)),
            rmse = sqrt(mean_or_na(e^2)),
            coverage = mean_or_na(data$covered),
            within_domain_rate = mean(data$within_domain %in% TRUE),
            failure_rate = mean(data$failed %in% TRUE),
            stringsAsFactors = FALSE
        )
    )
}))
rownames(by_design) <- NULL
write_tsv(by_design, file.path(cli$output_dir, "bslmm-validation-performance-by-design.tsv"))

## Stratified audits requested by the PI: null/low-PVE and high-dimensional.
high_dim_min <- as_num(settings$high_dim_num_snps_min %||% "5000",
                       "high_dim_num_snps_min")
audit_slice <- function(data, label) {
    avail <- is.finite(data$bslmm_pve) & !(data$failed %in% TRUE)
    e <- data$error
    data.frame(
        slice = label,
        n = nrow(data),
        n_estimable = sum(avail),
        bias = mean_or_na(e),
        absolute_mean_bias = abs(mean_or_na(e)),
        rmse = sqrt(mean_or_na(e^2)),
        null_mean_estimated_h2 = if (any(data$true_h2 == 0)) {
            mean_or_na(data$bslmm_pve[data$true_h2 == 0])
        } else {
            NA_real_
        },
        coverage = mean_or_na(data$covered),
        within_domain_rate = mean(data$within_domain %in% TRUE),
        failure_rate = mean(data$failed %in% TRUE),
        stringsAsFactors = FALSE
    )
}
audits <- rbind(
    audit_slice(evaluation, "overall"),
    audit_slice(evaluation[evaluation$true_h2 == 0, , drop = FALSE], "null_h2_0"),
    audit_slice(evaluation[evaluation$true_h2 <= low_pve_max, , drop = FALSE],
                paste0("low_pve_le_", low_pve_max)),
    audit_slice(evaluation[evaluation$num_snps >= high_dim_min, , drop = FALSE],
                paste0("high_dim_p_ge_", high_dim_min)),
    audit_slice(
        evaluation[evaluation$num_snps >= high_dim_min &
                       evaluation$architecture == "polygenic", , drop = FALSE],
        paste0("high_dim_polygenic_p_ge_", high_dim_min)
    ),
    audit_slice(
        evaluation[evaluation$num_snps >= high_dim_min & evaluation$true_h2 == 0,
                   , drop = FALSE],
        paste0("high_dim_null_p_ge_", high_dim_min)
    )
)
write_tsv(audits, file.path(cli$output_dir, "bslmm-validation-stratified-audits.tsv"))

## Module 02 acceptance criteria: hard biological-calibration gates plus the
## nonbinding aggregate-RMSE guardrail.
performance <- overall
if (!"gate_role" %in% names(criteria)) criteria$gate_role <- "hard"
if (!"gate_version" %in% names(criteria)) criteria$gate_version <- "legacy"
if (any(!criteria$gate_role %in% c("hard", "guardrail"))) {
    stop("gate_role must be hard or guardrail")
}
results <- lapply(seq_len(nrow(criteria)), function(i) {
    metric <- criteria$metric[[i]]
    if (!metric %in% names(performance)) {
        stop("Performance metric not found for BSLMM evaluation: ", metric)
    }
    observed <- as.numeric(performance[[metric]][[1L]])
    threshold <- as.numeric(criteria$threshold[[i]])
    comparison <- criteria$comparison[[i]]
    passed <- switch(
        comparison,
        less_than_or_equal = observed <= threshold,
        greater_than_or_equal = observed >= threshold,
        stop("Unknown comparison: ", comparison)
    )
    data.frame(
        gate_version = gate_version,
        metric = metric,
        observed = observed,
        comparison = comparison,
        threshold = threshold,
        gate_role = criteria$gate_role[[i]],
        passed = isTRUE(passed),
        rationale = criteria$rationale[[i]],
        stringsAsFactors = FALSE
    )
})
results <- do.call(rbind, results)

## Extra absolute gates analogous to observed-run rejection of EN.
min_within <- as_num(settings$min_within_domain_rate %||% "0.90",
                     "min_within_domain_rate")
max_fail <- as_num(settings$max_computational_failure_rate %||% "0",
                   "max_computational_failure_rate")
extra <- rbind(
    data.frame(
        gate_version = gate_version,
        metric = "within_domain_rate",
        observed = overall$within_domain_rate[[1L]],
        comparison = "greater_than_or_equal",
        threshold = min_within,
        gate_role = "hard",
        passed = isTRUE(overall$within_domain_rate[[1L]] >= min_within),
        rationale = "Analog of observed EN domain gate (max 10% outside)",
        stringsAsFactors = FALSE
    ),
    data.frame(
        gate_version = gate_version,
        metric = "computational_failure_rate",
        observed = overall$computational_failure_rate[[1L]],
        comparison = "less_than_or_equal",
        threshold = max_fail,
        gate_role = "hard",
        passed = isTRUE(overall$computational_failure_rate[[1L]] <= max_fail),
        rationale = "Zero unexplained computational failures",
        stringsAsFactors = FALSE
    ),
    data.frame(
        gate_version = gate_version,
        metric = "task_reconciliation_complete",
        observed = as.numeric(reconciliation$complete[[1L]]),
        comparison = "greater_than_or_equal",
        threshold = 1,
        gate_role = "hard",
        passed = isTRUE(reconciliation$complete[[1L]]),
        rationale = "Every expected validation scenario must complete",
        stringsAsFactors = FALSE
    )
)
results <- rbind(results, extra)
write_tsv(results, file.path(cli$output_dir, "bslmm-validation-acceptance-results.tsv"))

hard <- results$gate_role == "hard"
all_pass <- all(results$passed[hard])
decision <- data.frame(
    decision = if (all_pass) {
        "PASS_BSLMM_ELIGIBLE_TO_REPLACE_EN"
    } else {
        "FAIL_KEEP_EN_DO_NOT_REPLACE"
    },
    n_criteria = nrow(results),
    n_passed = sum(results$passed),
    n_hard_criteria = sum(hard),
    n_hard_passed = sum(results$passed[hard]),
    rmse = overall$rmse,
    rmse_guardrail_met = results$passed[results$metric == "rmse"],
    absolute_mean_bias = overall$absolute_mean_bias,
    null_mean_estimated_h2 = overall$null_mean_estimated_h2,
    low_h2_absolute_mean_bias = overall$low_h2_absolute_mean_bias,
    max_absolute_low_h2_level_bias = overall$max_absolute_low_h2_level_bias,
    max_absolute_h2_level_bias = overall$max_absolute_h2_level_bias,
    calibration_interval_coverage = overall$calibration_interval_coverage,
    null_type1_error = overall$null_type1_error,
    spearman_truth_estimate = overall$spearman_truth_estimate,
    within_domain_rate = overall$within_domain_rate,
    computational_failure_rate = overall$computational_failure_rate,
    note = paste(
        "Observed Module 02 production must not change unless decision is PASS.",
        "Recalibrated EN is optional and not required for this gate."
    ),
    stringsAsFactors = FALSE
)
write_tsv(decision, file.path(cli$output_dir, "bslmm-validation-decision.tsv"))

cat("\n=== BSLMM validation acceptance\n")
print(results, row.names = FALSE)
cat("\n=== decision\n")
print(decision, row.names = FALSE)
cat("\n=== stratified audits\n")
print(audits, row.names = FALSE)

if (as_bool(cli$fail_on_rejection, "fail_on_rejection") && !all_pass) {
    quit(save = "no", status = 2L)
}
