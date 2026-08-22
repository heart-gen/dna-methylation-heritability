#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))
source(file.path(dirname(script_path), "bslmm_pilot_functions.R"))

cli <- parse_cli(list(
    input_dir = "",
    manifest = "",
    config = "",
    output_dir = ""
))
if (!nzchar(cli$input_dir) || !nzchar(cli$manifest) ||
    !nzchar(cli$config) || !nzchar(cli$output_dir)) {
    stop("--input_dir, --manifest, --config, and --output_dir are required")
}
settings <- read_pilot_settings(cli$config)
manifest <- read_tsv(cli$manifest)
files <- list.files(
    cli$input_dir, pattern = "^scenario-[0-9]+\\.tsv$", full.names = TRUE
)
if (!length(files)) stop("No scenario result files found in ", cli$input_dir)
tables <- lapply(files, read_tsv)
columns <- Reduce(union, lapply(tables, names))
tables <- lapply(tables, function(d) {
    for (missing in setdiff(columns, names(d))) d[[missing]] <- NA
    d[, columns, drop = FALSE]
})
results <- do.call(rbind, tables)
results <- results[order(results$scenario_id), , drop = FALSE]
write_tsv(results, file.path(cli$output_dir, "bslmm-en-pilot-all-scenarios.tsv"))

expected <- nrow(manifest)
completed <- nrow(results)
missing_ids <- setdiff(manifest$scenario_id, results$scenario_id)
reconciliation <- data.frame(
    expected_scenarios = expected,
    completed_scenarios = completed,
    missing_scenarios = length(missing_ids),
    complete = length(missing_ids) == 0L,
    stringsAsFactors = FALSE
)
write_tsv(reconciliation, file.path(cli$output_dir, "task-reconciliation.tsv"))
if (length(missing_ids)) {
    write_tsv(
        data.frame(scenario_id = missing_ids),
        file.path(cli$output_dir, "missing-scenarios.tsv")
    )
}

metric_block <- function(estimate, truth, covered, failed, boundary = NULL) {
    keep <- is.finite(estimate) & is.finite(truth) & !failed
    err <- estimate[keep] - truth[keep]
    data.frame(
        n = length(truth),
        n_estimable = sum(keep),
        rmse = if (any(keep)) sqrt(mean(err^2)) else NA_real_,
        bias = if (any(keep)) mean(err) else NA_real_,
        absolute_mean_bias = if (any(keep)) abs(mean(err)) else NA_real_,
        coverage = if (length(covered)) mean(covered, na.rm = TRUE) else NA_real_,
        failure_rate = mean(failed, na.rm = TRUE),
        boundary_or_extrap_rate = if (is.null(boundary)) {
            NA_real_
        } else {
            mean(boundary, na.rm = TRUE)
        },
        stringsAsFactors = FALSE
    )
}

summarize_slice <- function(d, label) {
    en <- metric_block(
        d$h2_en_calibrated, d$true_h2, d$en_covered, d$en_failed,
        d$en_boundary_or_extrap
    )
    names(en) <- paste0("en_", names(en))
    bslmm <- metric_block(
        d$bslmm_pve, d$true_h2, d$bslmm_covered, d$bslmm_failed, NULL
    )
    names(bslmm) <- paste0("bslmm_", names(bslmm))
    cbind(data.frame(slice = label, stringsAsFactors = FALSE), en, bslmm)
}

overall <- summarize_slice(results, "overall")
by_p <- do.call(rbind, lapply(split(results, results$num_snps), function(d) {
    summarize_slice(d, paste0("num_snps=", d$num_snps[[1L]]))
}))
by_arch <- do.call(rbind, lapply(split(results, results$architecture), function(d) {
    summarize_slice(d, paste0("architecture=", d$architecture[[1L]]))
}))
by_p_arch <- do.call(rbind, lapply(
    split(results, list(results$num_snps, results$architecture), drop = TRUE),
    function(d) {
        summarize_slice(
            d,
            paste0(
                "num_snps=", d$num_snps[[1L]],
                "|architecture=", d$architecture[[1L]]
            )
        )
    }
))
summary <- rbind(overall, by_p, by_arch, by_p_arch)
write_tsv(summary, file.path(cli$output_dir, "bslmm-en-pilot-metrics.tsv"))

stress_arch <- settings$stress_architecture
stress_p_min <- as_num(settings$stress_num_snps_min, "stress_num_snps_min")
stress_h2_min <- as_num(settings$stress_h2_min, "stress_h2_min")
stress <- results[
    results$architecture == stress_arch &
        results$num_snps >= stress_p_min &
        results$true_h2 >= stress_h2_min,
    ,
    drop = FALSE
]
stress_metrics <- summarize_slice(stress, "stress_polygenic_high_p_high_h2")
write_tsv(stress_metrics, file.path(cli$output_dir, "bslmm-en-pilot-stress.tsv"))

rmse_rule <- as_num(
    settings$rmse_relative_improvement_min, "rmse_relative_improvement_min"
)
bias_rule <- as_num(
    settings$stress_bias_improvement_min, "stress_bias_improvement_min"
)
fail_delta <- as_num(
    settings$max_failure_rate_delta, "max_failure_rate_delta"
)

en_rmse <- overall$en_rmse[[1L]]
bslmm_rmse <- overall$bslmm_rmse[[1L]]
rmse_rel_improvement <- if (is.finite(en_rmse) && en_rmse > 0) {
    (en_rmse - bslmm_rmse) / en_rmse
} else {
    NA_real_
}
en_stress_bias <- stress_metrics$en_bias[[1L]]
bslmm_stress_bias <- stress_metrics$bslmm_bias[[1L]]
## Downward bias is negative. Improvement means moving toward zero from below,
## i.e. bslmm_bias - en_bias when both are negative, or reduction in |bias|.
stress_abs_bias_improvement <- if (
    is.finite(en_stress_bias) && is.finite(bslmm_stress_bias)
) {
    abs(en_stress_bias) - abs(bslmm_stress_bias)
} else {
    NA_real_
}
failure_delta <- overall$bslmm_failure_rate[[1L]] - overall$en_failure_rate[[1L]]

go_rmse <- is.finite(rmse_rel_improvement) && rmse_rel_improvement >= rmse_rule
go_bias <- is.finite(stress_abs_bias_improvement) &&
    stress_abs_bias_improvement >= bias_rule
go_fail <- is.finite(failure_delta) && failure_delta <= fail_delta
go <- isTRUE(reconciliation$complete[[1L]]) && go_rmse && go_bias && go_fail

decision <- data.frame(
    complete = reconciliation$complete[[1L]],
    en_rmse = en_rmse,
    bslmm_rmse = bslmm_rmse,
    rmse_relative_improvement = rmse_rel_improvement,
    rmse_rule = rmse_rule,
    rmse_pass = go_rmse,
    stress_n = nrow(stress),
    en_stress_bias = en_stress_bias,
    bslmm_stress_bias = bslmm_stress_bias,
    stress_abs_bias_improvement = stress_abs_bias_improvement,
    stress_bias_rule = bias_rule,
    stress_bias_pass = go_bias,
    en_failure_rate = overall$en_failure_rate[[1L]],
    bslmm_failure_rate = overall$bslmm_failure_rate[[1L]],
    failure_rate_delta = failure_delta,
    failure_delta_rule = fail_delta,
    failure_pass = go_fail,
    decision = if (go) "GO_EXPAND_VALIDATION_GRID" else "NO_GO_STOP",
    stringsAsFactors = FALSE
)
write_tsv(decision, file.path(cli$output_dir, "bslmm-en-pilot-decision.tsv"))

cat("\n=== paired pilot decision\n")
print(decision, row.names = FALSE)
if (!go) quit(save = "no", status = 2L)
