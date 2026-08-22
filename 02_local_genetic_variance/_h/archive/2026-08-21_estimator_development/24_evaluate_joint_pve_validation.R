#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))
source(file.path(dirname(script_path), "joint_pve_functions.R"))

cli <- parse_cli(list(
    input_dir = "", manifest = "", config = "", criteria = "",
    model = "", model_sha256 = "", output_dir = "",
    fail_on_rejection = "TRUE"
))
required <- c("input_dir", "manifest", "config", "criteria", "model",
              "model_sha256", "output_dir")
if (any(!nzchar(unlist(cli[required])))) stop("All required arguments must be set")
dir.create(cli$output_dir, recursive = TRUE, showWarnings = FALSE)
expected_sha <- trimws(readLines(cli$model_sha256, warn = FALSE)[[1L]])
observed_sha <- sub(" .*$", "", system2("sha256sum", cli$model, stdout = TRUE)[[1L]])
if (!identical(expected_sha, observed_sha)) stop("Joint model checksum mismatch")
model <- readRDS(cli$model)
settings <- read_joint_settings(cli$config)
criteria <- read_tsv(cli$criteria)
if (!identical(model$gate_version, unique(criteria$gate_version))) {
    stop("Model/acceptance gate version mismatch")
}
manifest <- read_tsv(cli$manifest)
files <- list.files(cli$input_dir, pattern = "^scenario-[0-9]+\\.tsv$",
                    full.names = TRUE)
data <- if (length(files)) do.call(rbind, lapply(files, read_tsv)) else
    data.frame()
duplicate_ids <- if (nrow(data)) sum(duplicated(data$scenario_id)) else 0L
missing_ids <- if (nrow(data)) setdiff(manifest$scenario_id, data$scenario_id) else
    manifest$scenario_id
extra_ids <- if (nrow(data)) setdiff(data$scenario_id, manifest$scenario_id) else
    integer()
if (duplicate_ids || length(missing_ids) || length(extra_ids)) {
    stop("Validation reconciliation failed before scientific evaluation")
}
data <- data[match(manifest$scenario_id, data$scenario_id), , drop = FALSE]
prediction <- predict_joint_pve(model, data)
evaluation <- cbind(data, prediction)
performance <- joint_validation_metrics(evaluation, expected_n = nrow(manifest))
acceptance <- evaluate_joint_criteria(performance$metrics, criteria)
hard <- acceptance$gate_role == "hard"
passed <- all(acceptance$passed[hard])

write_tsv(evaluation, file.path(cli$output_dir, "joint-pve-validation-estimates.tsv"))
write_tsv(performance$by_level,
          file.path(cli$output_dir, "joint-pve-validation-by-h2.tsv"))
write_tsv(acceptance,
          file.path(cli$output_dir, "joint-pve-validation-acceptance.tsv"))
overall <- data.frame(
    metric = names(performance$metrics),
    value = as.numeric(performance$metrics),
    stringsAsFactors = FALSE
)
write_tsv(overall, file.path(cli$output_dir, "joint-pve-validation-overall.tsv"))

strata <- c("architecture", "n", "num_snps", "ld_rho")
stratified <- do.call(rbind, lapply(strata, function(variable) {
    groups <- split(evaluation, evaluation[[variable]])
    do.call(rbind, lapply(names(groups), function(level) {
        d <- groups[[level]]
        err <- d$pve_cis_joint_calibrated - d$true_h2
        data.frame(
            variable = variable, level = level, n = nrow(d),
            bias = mean(err), rmse = sqrt(mean(err^2)),
            spearman = suppressWarnings(stats::cor(
                d$true_h2, d$pve_cis_joint_calibrated, method = "spearman"
            )),
            lower_boundary_rate = mean(d$pve_lower_boundary_hit),
            upper_boundary_rate = mean(d$pve_upper_boundary_hit),
            stringsAsFactors = FALSE
        )
    }))
}))
write_tsv(stratified,
          file.path(cli$output_dir, "joint-pve-validation-stratified.tsv"))

decision <- data.frame(
    experiment = "final_joint_pve_v1",
    decision = if (passed) "PASS_PROMOTE_ABSOLUTE_PVE" else
        "FAIL_PIVOT_TO_RELATIVE_LOCAL_GENETIC_CONTROL",
    hard_criteria_passed = sum(acceptance$passed[hard]),
    hard_criteria_total = sum(hard),
    model_sha256 = observed_sha,
    validation_scenarios = nrow(manifest),
    stringsAsFactors = FALSE
)
write_tsv(decision, file.path(cli$output_dir, "joint-pve-terminal-decision.tsv"))

failed_metrics <- acceptance$metric[hard & !acceptance$passed]
decision_md <- c(
    "# Module 02 terminal absolute-PVE decision",
    "",
    paste0("**Decision:** `", decision$decision, "`"),
    "",
    paste0("Independent validation scenarios: ", nrow(manifest), "."),
    paste0("Hard criteria passed: ", sum(acceptance$passed[hard]), "/",
           sum(hard), "."),
    if (passed) {
        paste(
            "The locked joint estimator may be described as simulation-calibrated",
            "local cis-SNP PVE. Observed-data promotion still requires the accepted",
            "VMR catalog and a separately reconciled observed run."
        )
    } else {
        paste(
            "Absolute locus-level PVE is not identifiable with sufficient reliability",
            "at these sample sizes. Estimator development ends; Module 02 and the",
            "manuscript move to a continuous relative/local-genetic-control axis."
        )
    },
    if (length(failed_metrics)) paste0("Failed hard metrics: ",
        paste(failed_metrics, collapse = ", "), ".") else
        "Failed hard metrics: none.",
    "",
    paste0("Frozen model SHA-256: `", observed_sha, "`.")
)
writeLines(decision_md, file.path(cli$output_dir, "MODULE02_TERMINAL_DECISION.md"))
cat("Final joint PVE decision:", decision$decision, "\n")
if (!passed && as_bool(cli$fail_on_rejection, "fail_on_rejection")) quit(status = 2L)
