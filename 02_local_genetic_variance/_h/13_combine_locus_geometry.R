#!/usr/bin/env Rscript

## Locus-geometry scan, stage C: combine and reconcile the task universe.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(run_dir = ""))
if (!nzchar(cli$run_dir)) stop("--run-dir is required")
run_dir <- normalizePath(cli$run_dir)
tasks <- read_tsv(file.path(run_dir, "config", "task-manifest.tsv"))
row_dir <- file.path(run_dir, "results", "task_rows")
row_files <- list.files(row_dir, pattern = "^vmr-[0-9]{7}\\.tsv$",
                        full.names = TRUE)
if (!length(row_files)) stop("No geometry rows to combine")
rows <- do.call(rbind, lapply(row_files, read_tsv))
rows <- rows[order(rows$task_id), , drop = FALSE]
combined_dir <- file.path(run_dir, "results", "combined")
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(rows, file.path(combined_dir, "locus-geometry.tsv"))

expected <- as.integer(tasks$task_id)
seen <- as.integer(rows$task_id)
recon <- data.frame(
    expected = length(expected), task_files = length(row_files),
    unique_task_rows = length(unique(seen)),
    completed = sum(rows$terminal_status %in% "completed"),
    qc_failed = sum(rows$terminal_status %in% "qc_failed"),
    excluded = sum(rows$terminal_status %in% "excluded"),
    computational_failure = sum(rows$terminal_status %in% "computational_failure"),
    duplicate_task_ids = sum(duplicated(seen)),
    unexpected_task_ids = length(setdiff(seen, expected)),
    unaccounted = length(setdiff(expected, seen)),
    stringsAsFactors = FALSE
)
write_tsv(recon, file.path(combined_dir, "geometry-reconciliation.tsv"))
missing <- setdiff(expected, seen)
write_tsv(data.frame(task_id = missing),
          file.path(combined_dir, "missing-task-ids.tsv"))
cat("Combined", nrow(rows), "geometry rows;", recon$completed, "complete;",
    length(missing), "unaccounted\n")
if (length(missing) || recon$duplicate_task_ids ||
    recon$unexpected_task_ids) {
    quit(save = "no", status = 1L)
}
