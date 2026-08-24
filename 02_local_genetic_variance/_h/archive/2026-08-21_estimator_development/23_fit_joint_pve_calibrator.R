#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))
source(file.path(dirname(script_path), "joint_pve_functions.R"))

cli <- parse_cli(list(
    input_dir = "", manifest = "", config = "", output_dir = ""
))
if (any(!nzchar(unlist(cli)))) stop("All arguments are required")
dir.create(cli$output_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- read_tsv(cli$manifest)
settings <- read_joint_settings(cli$config)
files <- list.files(cli$input_dir, pattern = "^scenario-[0-9]+\\.tsv$",
                    full.names = TRUE)
if (!length(files)) stop("No development feature files found")
data <- do.call(rbind, lapply(files, read_tsv))
if (anyDuplicated(data$scenario_id)) stop("Duplicate development scenario outputs")
missing_ids <- setdiff(manifest$scenario_id, data$scenario_id)
extra_ids <- setdiff(data$scenario_id, manifest$scenario_id)
if (length(missing_ids) || length(extra_ids) || nrow(data) != nrow(manifest)) {
    stop("Development reconciliation failed: expected=", nrow(manifest),
         " completed=", nrow(data), " missing=", length(missing_ids),
         " extra=", length(extra_ids))
}
data <- data[match(manifest$scenario_id, data$scenario_id), , drop = FALSE]
fit_reps <- as.integer(split_numeric(settings$development_fit_replicates))
cal_reps <- as.integer(split_numeric(settings$development_calibration_replicates))
if (length(intersect(fit_reps, cal_reps))) stop("Fit/calibration replicates overlap")
fit_data <- data[data$replicate %in% fit_reps, , drop = FALSE]
cal_data <- data[data$replicate %in% cal_reps, , drop = FALSE]
if (nrow(fit_data) + nrow(cal_data) != nrow(data)) {
    stop("Some development replicates are not assigned")
}
model <- fit_joint_pve_model(fit_data, cal_data, settings)
model$training_manifest_rows <- nrow(manifest)
model$development_source <- settings$development_bslmm_run_relpath
model$fit_created_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)

model_path <- file.path(cli$output_dir, "joint-pve-calibrator.rds")
saveRDS(model, model_path, version = 3)
sha <- system2("sha256sum", model_path, stdout = TRUE)
sha <- sub(" .*$", "", sha[[1L]])
writeLines(sha, file.path(cli$output_dir, "joint-pve-calibrator.sha256"))

development_prediction <- rbind(
    cbind(subset = "fit", fit_data, predict_joint_pve(model, fit_data)),
    cbind(subset = "calibration", cal_data,
          predict_joint_pve(model, cal_data))
)
write_tsv(data, file.path(cli$output_dir, "development-features.tsv"))
write_tsv(development_prediction,
          file.path(cli$output_dir, "development-predictions.tsv"))

coef_matrix <- as.matrix(stats::coef(model$glmnet_fit, s = model$lambda))
coefficients <- data.frame(
    term = rownames(coef_matrix),
    coefficient = as.numeric(coef_matrix[, 1L]),
    group = ifelse(
        startsWith(rownames(coef_matrix), "signal__"),
        sub("^signal__(.*?)__.*$", "\\1", rownames(coef_matrix), perl = TRUE),
        ifelse(startsWith(rownames(coef_matrix), "design__"), "design", "intercept")
    ),
    stringsAsFactors = FALSE
)
write_tsv(coefficients, file.path(cli$output_dir, "joint-pve-coefficients.tsv"))
metadata <- data.frame(
    field = c(
        "family", "gate_version", "ridge_lambda", "fit_rows",
        "calibration_rows", "null_calibration_rows", "conformal_q",
        "null_cutoff", "model_sha256"
    ),
    value = c(
        model$family, model$gate_version, model$lambda,
        model$development_counts[["fit"]],
        model$development_counts[["calibration"]],
        model$development_counts[["null"]], model$conformal_q,
        model$null_cutoff, sha
    ),
    stringsAsFactors = FALSE
)
write_tsv(metadata, file.path(cli$output_dir, "joint-pve-model-metadata.tsv"))
cat("Frozen joint PVE model:", model_path, "\nSHA-256:", sha, "\n")
