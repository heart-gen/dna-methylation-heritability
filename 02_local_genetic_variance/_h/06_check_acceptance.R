#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    performance = file.path(dirname(script_path), "..", "_m", "evaluation", "calibration-performance-overall.tsv"),
    criteria = file.path(dirname(script_path), "..", "config", "acceptance-criteria.tsv"),
    output = file.path(dirname(script_path), "..", "_m", "evaluation", "acceptance-results.tsv"),
    fail_on_rejection = "FALSE"
))
performance <- read_tsv(cli$performance)
criteria <- read_tsv(cli$criteria)
if (nrow(performance) != 1L) stop("Overall performance file must contain one row")

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
        metric = metric,
        observed = observed,
        comparison = comparison,
        threshold = threshold,
        passed = isTRUE(passed),
        rationale = criteria$rationale[[i]],
        stringsAsFactors = FALSE
    )
})
results <- do.call(rbind, results)
write_tsv(results, cli$output)
cat("Acceptance criteria passed:", sum(results$passed), "of", nrow(results), "\n")
if (as_bool(cli$fail_on_rejection, "fail_on_rejection") && !all(results$passed)) {
    quit(save = "no", status = 2L)
}
