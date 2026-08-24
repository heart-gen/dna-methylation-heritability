#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    performance = file.path(dirname(script_path), "..", "_m", "evaluation", "calibration-performance-overall.tsv"),
    criteria = file.path(dirname(script_path), "..", "config", "acceptance-criteria.tsv"),
    output = file.path(dirname(script_path), "..", "_m", "evaluation", "acceptance-results.tsv"),
    model = "",
    fail_on_rejection = "FALSE"
))
performance <- read_tsv(cli$performance)
criteria <- read_tsv(cli$criteria)
if (nrow(performance) != 1L) stop("Overall performance file must contain one row")
if (!"gate_role" %in% names(criteria)) criteria$gate_role <- "hard"
if (!"gate_version" %in% names(criteria)) criteria$gate_version <- "legacy"
if (any(!criteria$gate_role %in% c("hard", "guardrail"))) {
    stop("gate_role must be hard or guardrail")
}
gate_versions <- unique(criteria$gate_version)
if (length(gate_versions) != 1L || !nzchar(gate_versions[[1L]])) {
    stop("Acceptance criteria must declare exactly one non-empty gate_version")
}
gate_version <- gate_versions[[1L]]

results <- lapply(seq_len(nrow(criteria)), function(i) {
    metric <- criteria$metric[[i]]
    if (!metric %in% names(performance)) stop("Performance metric not found: ", metric)
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

# Passing the independent numerical criteria is necessary but not sufficient:
# the model-fitting stage also requires a weight satisfying every locked
# internal-development constraint. Production callers pass --model so that a
# fallback model retained for diagnostics cannot be applied to observed data.
if (nzchar(cli$model)) {
    if (!file.exists(cli$model)) stop("Calibration model not found: ", cli$model)
    model <- readRDS(cli$model)
    model_gate_version <- model$acceptance_gate_version
    if (is.null(model_gate_version) || length(model_gate_version) != 1L) {
        model_gate_version <- "missing"
    }
    internal_status <- model$internal_selection_status
    if (is.null(internal_status) || length(internal_status) != 1L) {
        internal_status <- "missing"
    }
    internal_passed <- identical(internal_status, "constraints_satisfied") &&
        identical(model_gate_version, gate_version)
    results <- rbind(
        results,
        data.frame(
            gate_version = gate_version,
            metric = "internal_development_constraints",
            observed = as.numeric(internal_passed),
            comparison = "greater_than_or_equal",
            threshold = 1,
            gate_role = "hard",
            passed = internal_passed,
            rationale = paste0(
                "The calibration model must record constraints_satisfied; ",
                "recorded status was ", internal_status,
                ". Model gate version must be ", gate_version,
                "; recorded version was ", model_gate_version
            ),
            stringsAsFactors = FALSE
        )
    )
}
write_tsv(results, cli$output)
hard <- results$gate_role == "hard"
hard_pass <- all(results$passed[hard])
cat("Hard acceptance criteria passed:", sum(results$passed[hard]), "of", sum(hard), "\n")
if (any(!hard)) {
    cat("Guardrails met:", sum(results$passed[!hard]), "of", sum(!hard), "\n")
}
if (as_bool(cli$fail_on_rejection, "fail_on_rejection") && !hard_pass) {
    quit(save = "no", status = 2L)
}
