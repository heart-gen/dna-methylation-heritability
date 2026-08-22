#!/usr/bin/env Rscript

## Stage 02: reconcile the Stage 01 task universe and create one combined row
## per expected VMR. Missing scheduler outputs become explicit computational
## failures rather than silently disappearing from denominators.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(run_dir = ""))
if (!nzchar(cli$run_dir)) stop("--run-dir is required")
run_dir <- normalizePath(cli$run_dir)
tasks <- read_tsv(file.path(run_dir, "config", "task-manifest.tsv"))
files <- list.files(file.path(run_dir, "results", "task_rows"),
                    pattern = "^vmr-[0-9]{7}\\.tsv$", full.names = TRUE)
rows <- if (length(files)) do.call(rbind, lapply(files, read_tsv)) else NULL

duplicate_ids <- if (is.null(rows)) integer() else
    unique(rows$task_id[duplicated(rows$task_id) |
                        duplicated(rows$task_id, fromLast = TRUE)])
unexpected_ids <- if (is.null(rows)) integer() else
    setdiff(rows$task_id, tasks$task_id)
observed_ids <- if (is.null(rows)) integer() else unique(rows$task_id)
missing_ids <- setdiff(tasks$task_id, observed_ids)

reconciliation <- data.frame(
    expected = nrow(tasks),
    task_files = length(files),
    unique_task_rows = length(observed_ids),
    completed = if (is.null(rows)) 0L else sum(rows$terminal_status == "completed"),
    excluded = if (is.null(rows)) 0L else sum(rows$terminal_status == "excluded"),
    qc_failed = if (is.null(rows)) 0L else sum(rows$terminal_status == "qc_failed"),
    computational_failure = if (is.null(rows)) 0L else
        sum(rows$terminal_status == "computational_failure"),
    unaccounted = length(missing_ids),
    duplicate_task_ids = length(duplicate_ids),
    unexpected_task_ids = length(unexpected_ids),
    stringsAsFactors = FALSE
)
combined_dir <- file.path(run_dir, "results", "combined")
write_tsv(reconciliation,
          file.path(combined_dir, "task-reconciliation.tsv"))
if (length(missing_ids)) {
    write_tsv(data.frame(task_id = missing_ids),
              file.path(combined_dir, "missing-task-ids.tsv"))
}
if (length(duplicate_ids)) {
    write_tsv(data.frame(task_id = duplicate_ids),
              file.path(combined_dir, "duplicate-task-ids.tsv"))
}
if (length(unexpected_ids)) {
    write_tsv(data.frame(task_id = unexpected_ids),
              file.path(combined_dir, "unexpected-task-ids.tsv"))
}
if (length(duplicate_ids) || length(unexpected_ids)) {
    stop("Task reconciliation found duplicate or unexpected task IDs")
}
if (is.null(rows)) stop("No Stage 01 task rows were produced")

## Materialize missing rows with the same schema so every expected VMR remains
## visible downstream. Stage 05 will fail acceptance when any are present.
if (length(missing_ids)) {
    manifest <- read_tsv(file.path(run_dir, "manifest.tsv"))
    upstream_vmr_run_id <- manifest$value[
        manifest$field == "upstream_vmr_run_id"
    ]
    if (length(upstream_vmr_run_id) != 1L) {
        stop("Manifest lacks one upstream_vmr_run_id")
    }
    missing <- tasks[tasks$task_id %in% missing_ids, , drop = FALSE]
    template <- rows[rep(1L, nrow(missing)), , drop = FALSE]
    template[,] <- NA
    shared <- intersect(names(template), names(missing))
    template[shared] <- missing[shared]
    template$population <- missing$cohort
    template$upstream_vmr_run_id <- upstream_vmr_run_id
    template$feature_complete <- FALSE
    template$computational_failure <- TRUE
    template$terminal_status <- "computational_failure"
    template$feature_error <- "missing Stage 01 task output"
    rows <- rbind(rows, template)
}
rows <- rows[order(rows$task_id), , drop = FALSE]
if (!identical(as.integer(rows$task_id), as.integer(tasks$task_id))) {
    stop("Combined rows do not reconcile exactly to the ordered task manifest")
}
write_tsv(rows, file.path(combined_dir, "observed-joint-features.tsv"))
cat("Reconciled", nrow(rows), "observed VMR tasks\n")
