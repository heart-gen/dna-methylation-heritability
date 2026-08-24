#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    base_config = "",
    recovery_lock = "",
    output = ""
))
for (key in names(cli)) {
    if (!nzchar(cli[[key]])) stop("--", key, " is required")
}

base <- read_tsv(cli$base_config)
lock <- read_tsv(cli$recovery_lock)
if (!all(c("setting", "value", "description") %in% names(base))) {
    stop("Base config must contain setting, value, and description")
}
if (!all(c("field", "value") %in% names(lock))) {
    stop("Recovery lock must contain field and value")
}
lock_values <- setNames(as.character(lock$value), lock$field)
required <- c("validation_seed_offset", "validation_run_id",
              "allowed_base_config_overrides", "model_changed",
              "criteria_changed")
missing <- setdiff(required, names(lock_values))
if (length(missing)) stop("Recovery lock is missing: ", paste(missing, collapse = ", "))
if (!identical(lock_values[["model_changed"]], "FALSE") ||
    !identical(lock_values[["criteria_changed"]], "FALSE")) {
    stop("Recovery may not change the frozen model or criteria")
}
allowed <- trimws(strsplit(lock_values[["allowed_base_config_overrides"]],
                           ",", fixed = TRUE)[[1L]])
expected <- c("validation_seed_offset", "validation_run_id")
if (!setequal(allowed, expected)) {
    stop("Only validation_seed_offset and validation_run_id may change")
}
if (!all(expected %in% base$setting)) {
    stop("Base config lacks the allowed recovery fields")
}
original <- base
for (key in expected) {
    base$value[base$setting == key] <- lock_values[[key]]
}
changed <- base$setting[base$value != original$value]
if (!setequal(changed, expected)) {
    stop("Effective recovery config changed unexpected fields: ",
         paste(changed, collapse = ", "))
}
write_tsv(base, cli$output)
cat("Wrote recovery config with exactly these changes:",
    paste(sort(changed), collapse = ", "), "\n")
