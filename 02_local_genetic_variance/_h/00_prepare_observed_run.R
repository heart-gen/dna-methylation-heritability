#!/usr/bin/env Rscript

## Stage 00: open one immutable observed local-genetic-control run.
## This stage performs no estimation. It resolves and freezes the accepted VMR
## catalog, final joint model, task universe, settings, and checksums.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
h_dir <- dirname(script_path)
source(file.path(h_dir, "00_functions.R"))

repo_root <- normalizePath(file.path(h_dir, "..", ".."))
module_root <- file.path(repo_root, "02_local_genetic_variance")
cli <- parse_cli(list(
    run_id = "",
    cohort = "",
    region = "",
    vmr_run_id = "",
    model_run_id = "lgv-joint-pve-train-20260820",
    model_sha256 =
        "9f26c3273746fda85d9bbf21e224857db9a1ad79a521582a12f241854c03223a",
    smoke_n = "0",
    vmrs_per_chunk = "5",
    runs_root = ""
))

required_cli <- c("run_id", "cohort", "region", "vmr_run_id")
missing_cli <- required_cli[!vapply(cli[required_cli], nzchar, logical(1L))]
if (length(missing_cli)) {
    stop("Required argument(s): ", paste(missing_cli, collapse = ", "))
}
cohort <- cli$cohort
region <- tolower(cli$region)
if (!cohort %in% c("AA", "all_individuals")) stop("Unsupported cohort: ", cohort)
if (!region %in% c("caudate", "dlpfc", "hippocampus")) {
    stop("Unsupported region: ", region)
}
expected_run_pattern <- paste0(
    "^lgv-", cohort, "-", region, "-[0-9]{8}[a-z]?$"
)
if (!grepl(expected_run_pattern, cli$run_id)) {
    stop("run_id must match ", expected_run_pattern)
}

runs_root <- if (nzchar(cli$runs_root)) cli$runs_root else
    file.path(module_root, "_m", "runs")
dir.create(runs_root, recursive = TRUE, showWarnings = FALSE)
run_dir <- file.path(runs_root, cli$run_id)
if (file.exists(run_dir)) stop("Run directory already exists: ", run_dir)

vmr_run_dir <- file.path(
    repo_root, "01_vmr_catalog", "_m", "runs", cli$vmr_run_id
)
vmr_manifest_path <- file.path(vmr_run_dir, "manifest.tsv")
vmr_catalog_path <- file.path(vmr_run_dir, "vmr", "vmr_catalog.tsv")
for (path in c(vmr_manifest_path, vmr_catalog_path)) {
    if (!file.exists(path)) stop("Required Module 01 input is missing: ", path)
}

## Module 01 predates the shared `Accepted runs` parser. Its README is still
## the record of acceptance, so require the exact run row and its locked gate.
vmr_readme <- readLines(file.path(repo_root, "01_vmr_catalog", "README.md"),
                        warn = FALSE)
acceptance_row <- vmr_readme[grepl(
    paste0("| `", cli$vmr_run_id, "` |"), vmr_readme, fixed = TRUE
)]
if (length(acceptance_row) != 1L ||
    !grepl("all five pass", acceptance_row, fixed = TRUE)) {
    stop("Module 01 run is not recorded as passing all five gates: ",
         cli$vmr_run_id)
}

vmr_manifest <- read_tsv(vmr_manifest_path)
manifest_value <- function(field) {
    value <- vmr_manifest$value[vmr_manifest$field == field]
    if (length(value) != 1L) stop("Module 01 manifest lacks unique field: ", field)
    as.character(value[[1L]])
}
if (!identical(manifest_value("cohort"), cohort) ||
    !identical(tolower(manifest_value("region")), region)) {
    stop("Requested cohort/region does not match Module 01 manifest")
}
if (!identical(toupper(manifest_value("smoke_run")), "FALSE")) {
    stop("Observed production must start from a non-smoke Module 01 run")
}

model_dir <- file.path(
    module_root, "_m", "runs", cli$model_run_id, "combined"
)
model_path <- file.path(model_dir, "joint-pve-calibrator.rds")
development_features <- file.path(model_dir, "development-features.tsv")
for (path in c(model_path, development_features)) {
    if (!file.exists(path)) stop("Frozen joint-model input is missing: ", path)
}
sha256_file <- function(path) {
    out <- system2("sha256sum", normalizePath(path), stdout = TRUE)
    if (!length(out)) stop("Could not checksum ", path)
    tolower(sub(" .*$", "", out[[1L]]))
}
observed_model_sha <- sha256_file(model_path)
if (!identical(observed_model_sha, tolower(cli$model_sha256))) {
    stop("Frozen joint-model checksum mismatch\n  expected ", cli$model_sha256,
         "\n  observed ", observed_model_sha)
}

tasks <- read_tsv(vmr_catalog_path)
required_task <- c("chr", "start", "end", "n", "vmr_id", "vmr_set_id")
missing_task <- setdiff(required_task, names(tasks))
if (length(missing_task)) {
    stop("VMR catalog lacks: ", paste(missing_task, collapse = ", "))
}
if (anyDuplicated(tasks$vmr_id)) stop("Duplicate vmr_id in VMR catalog")
if (length(unique(tasks$vmr_set_id)) != 1L ||
    !identical(unique(tasks$vmr_set_id), manifest_value("vmr_set_id"))) {
    stop("VMR catalog and manifest vmr_set_id differ")
}
tasks <- data.frame(
    task_id = seq_len(nrow(tasks)),
    cohort = cohort,
    region = region,
    chrom = as.character(tasks$chr),
    start = as.integer(tasks$start),
    end = as.integer(tasks$end),
    n_cpgs = as.integer(tasks$n),
    vmr_id = as.character(tasks$vmr_id),
    vmr_set_id = as.character(tasks$vmr_set_id),
    stringsAsFactors = FALSE
)
smoke_n <- as_int(cli$smoke_n, "smoke_n")
if (smoke_n < 0L) stop("smoke_n cannot be negative")
smoke_run <- smoke_n > 0L
if (smoke_run) tasks <- tasks[seq_len(min(smoke_n, nrow(tasks))), , drop = FALSE]
tasks$task_id <- seq_len(nrow(tasks))
vmrs_per_chunk <- as_int(cli$vmrs_per_chunk, "vmrs_per_chunk")
if (vmrs_per_chunk < 1L) stop("vmrs_per_chunk must be positive")
chunk_manifest <- data.frame(
    chunk_id = ceiling(tasks$task_id / vmrs_per_chunk),
    task_id = tasks$task_id,
    vmr_id = tasks$vmr_id,
    stringsAsFactors = FALSE
)

joint_config <- file.path(module_root, "config", "joint-pve-20260820.tsv")
relative_config <- file.path(repo_root, "config", "local_genetic_control.yml")
threshold_config <- file.path(repo_root, "config", "thresholds.yml")
for (path in c(joint_config, relative_config, threshold_config)) {
    if (!file.exists(path)) stop("Locked configuration is missing: ", path)
}

dir.create(run_dir, recursive = TRUE)
for (subdir in c("config", "logs", "work", "results/task_rows",
                 "results/combined")) {
    dir.create(file.path(run_dir, subdir), recursive = TRUE)
}
write_tsv(tasks, file.path(run_dir, "config", "task-manifest.tsv"))
write_tsv(chunk_manifest, file.path(run_dir, "config", "chunk-manifest.tsv"))
invisible(file.copy(c(joint_config, relative_config, threshold_config),
                    file.path(run_dir, "config"), overwrite = FALSE))

git_commit <- tryCatch(
    system2("git", c("-C", repo_root, "rev-parse", "HEAD"), stdout = TRUE)[[1L]],
    error = function(e) NA_character_
)
manifest <- data.frame(
    field = c(
        "run_id", "analysis", "cohort", "region", "started_at",
        "git_commit", "smoke_run", "upstream_vmr_run_id", "vmr_set_id",
        "ordered_donor_checksum", "n_donors", "n_expected_tasks",
        "vmrs_per_chunk", "n_expected_chunks",
        "joint_model_run_id", "joint_model_path", "joint_model_sha256",
        "development_features_path", "config_joint_sha256",
        "config_relative_score_sha256", "config_thresholds_sha256"
    ),
    value = c(
        cli$run_id, "02_local_genetic_variance", cohort, region,
        format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), git_commit,
        toupper(as.character(smoke_run)), cli$vmr_run_id,
        manifest_value("vmr_set_id"), manifest_value("donor_checksum"),
        manifest_value("n_donors"), nrow(tasks), vmrs_per_chunk,
        length(unique(chunk_manifest$chunk_id)), cli$model_run_id,
        normalizePath(model_path), observed_model_sha,
        normalizePath(development_features), sha256_file(joint_config),
        sha256_file(relative_config), sha256_file(threshold_config)
    ),
    stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(run_dir, "manifest.tsv"))
cat(normalizePath(run_dir), "\n", sep = "")
