#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    initial_qc = "",
    replacement_qc = "",
    initial_estimates = "",
    replacement_estimates = "",
    output_dir = file.path(dirname(script_path), "..", "_m", "comparison")
))
for (key in c("initial_qc", "replacement_qc", "initial_estimates",
              "replacement_estimates")) {
    if (!nzchar(cli[[key]])) stop("--", gsub("_", "-", key), " is required")
}

initial_qc <- read_tsv(cli$initial_qc)
replacement_qc <- read_tsv(cli$replacement_qc)
initial <- read_tsv(cli$initial_estimates)
replacement <- read_tsv(cli$replacement_estimates)
banned <- c("r_squared_cv", "h2_unscaled")
if (length(intersect(banned, c(names(initial), names(replacement))))) {
    stop("A comparison input contains a banned legacy metric")
}

cell_key <- c("population", "region", "vmr_set_id")
for (data in list(initial_qc, replacement_qc)) {
    missing <- setdiff(c(cell_key, "within_domain_rate", "overall_qc_pass"),
                       names(data))
    if (length(missing)) stop("QC input lacks: ", paste(missing, collapse = ", "))
}
gate <- merge(
    initial_qc[, c(cell_key, "analyzed_tasks", "within_domain_rate",
                   "overall_qc_pass")],
    replacement_qc[, c(cell_key, "analyzed_tasks", "within_domain_rate",
                       "overall_qc_pass")],
    by = cell_key, suffixes = c("_initial", "_replacement"), all = TRUE
)
if (nrow(gate) != 6L || anyNA(gate$within_domain_rate_initial) ||
    anyNA(gate$within_domain_rate_replacement)) {
    stop("Initial and replacement QC tables do not reconcile to six cells")
}
gate$within_domain_rate_change <- gate$within_domain_rate_replacement -
    gate$within_domain_rate_initial
gate$gate_improved <- gate$within_domain_rate_change > 0
gate$replacement_gate_pass <- gate$within_domain_rate_replacement >= 0.90 &
    gate$overall_qc_pass_replacement

key <- c("population", "region", "vmr_set_id", "vmr_id")
matched <- merge(
    initial[, c(key, "rho2_oof", "h2_en_calibrated", "calibration_status")],
    replacement[, c(key, "rho2_oof", "h2_en_calibrated", "calibration_status")],
    by = key, suffixes = c("_initial", "_replacement")
)
if (!nrow(matched)) stop("No matched VMR estimates between runs")
matched$rho2_change <- matched$rho2_oof_replacement - matched$rho2_oof_initial
matched$h2_change <- matched$h2_en_calibrated_replacement -
    matched$h2_en_calibrated_initial
match_key <- interaction(matched$population, matched$region, drop = TRUE,
                         lex.order = TRUE)
matched_summary <- do.call(rbind, lapply(split(matched, match_key), function(x) {
    data.frame(
        population = x$population[[1L]],
        region = x$region[[1L]],
        matched_vmrs = nrow(x),
        rho2_spearman = suppressWarnings(stats::cor(
            x$rho2_oof_initial, x$rho2_oof_replacement,
            method = "spearman", use = "complete.obs"
        )),
        median_rho2_change = median(x$rho2_change, na.rm = TRUE),
        median_h2_change = median(x$h2_change, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
}))
rownames(matched_summary) <- NULL

decision <- data.frame(
    all_six_replacement_gates_pass = all(gate$replacement_gate_pass),
    all_six_within_domain_rates_improved = all(gate$gate_improved),
    remove_initial_archive = all(gate$replacement_gate_pass) &&
        all(gate$gate_improved),
    decision_rule = paste(
        "Remove only if all six replacement runs pass complete/failure/domain QC",
        "and each within-domain rate exceeds its initial value"
    ),
    stringsAsFactors = FALSE
)

dir.create(cli$output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(gate, file.path(cli$output_dir, "within-domain-gate-comparison.tsv"))
write_tsv(matched_summary,
          file.path(cli$output_dir, "matched-estimator-comparison-by-cell.tsv"))
write_tsv(decision, file.path(cli$output_dir, "archive-retirement-decision.tsv"))
cat("Archive retirement authorized:", decision$remove_initial_archive, "\n")
