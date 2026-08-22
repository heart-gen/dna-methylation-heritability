#!/usr/bin/env Rscript

## Locus-geometry scan, stage A: open a run and freeze its task universe.
##
## The observed-regime grid stratifies real loci on num_snps and p_eff. For the
## caudate AA cell those came from a completed production run, but the other
## five cohort-by-region cells have none, and running a full production pass
## purely to enable a diagnostic would be backwards. p_eff and ld_metric are
## pure genotype geometry: no phenotype, no elastic net, no BSLMM. This scan
## computes them directly from the accepted Module 01 catalog.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
h_dir <- dirname(script_path)
source(file.path(h_dir, "00_functions.R"))

repo_root <- normalizePath(file.path(h_dir, "..", ".."))
module_root <- file.path(repo_root, "02_local_genetic_variance")
cli <- parse_cli(list(
    run_id = "", cohort = "", region = "", vmr_run_id = "",
    vmrs_per_chunk = "25", runs_root = ""
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
if (!grepl(paste0("^lgv-geometry-", cohort, "-", region, "-[0-9]{8}[a-z]?$"),
           cli$run_id)) {
    stop("run_id must match ^lgv-geometry-", cohort, "-", region,
         "-[0-9]{8}[a-z]?$")
}

runs_root <- if (nzchar(cli$runs_root)) cli$runs_root else
    file.path(module_root, "_m", "runs")
run_dir <- file.path(runs_root, cli$run_id)
if (file.exists(run_dir)) stop("Run directory already exists: ", run_dir)

vmr_run_dir <- file.path(repo_root, "01_vmr_catalog", "_m", "runs",
                         cli$vmr_run_id)
vmr_manifest_path <- file.path(vmr_run_dir, "manifest.tsv")
vmr_catalog_path <- file.path(vmr_run_dir, "vmr", "vmr_catalog.tsv")
for (path in c(vmr_manifest_path, vmr_catalog_path)) {
    if (!file.exists(path)) stop("Required Module 01 input is missing: ", path)
}
vmr_readme <- readLines(file.path(repo_root, "01_vmr_catalog", "README.md"),
                        warn = FALSE)
acceptance_row <- vmr_readme[grepl(paste0("| `", cli$vmr_run_id, "` |"),
                                   vmr_readme, fixed = TRUE)]
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
    stop("Module 01 run does not match the requested cohort and region")
}

catalog <- read_tsv(vmr_catalog_path)
required_task <- c("chr", "start", "end", "n", "vmr_id", "vmr_set_id")
missing_task <- setdiff(required_task, names(catalog))
if (length(missing_task)) {
    stop("VMR catalog lacks: ", paste(missing_task, collapse = ", "))
}
if (anyDuplicated(catalog$vmr_id)) stop("Duplicate vmr_id in VMR catalog")
if (!identical(unique(catalog$vmr_set_id), manifest_value("vmr_set_id"))) {
    stop("VMR catalog and manifest vmr_set_id differ")
}
tasks <- data.frame(
    task_id = seq_len(nrow(catalog)),
    cohort = cohort, region = region,
    chrom = as.character(catalog$chr),
    start = as.integer(catalog$start),
    end = as.integer(catalog$end),
    n_cpgs = as.integer(catalog$n),
    vmr_id = as.character(catalog$vmr_id),
    vmr_set_id = as.character(catalog$vmr_set_id),
    stringsAsFactors = FALSE
)
vmrs_per_chunk <- as_int(cli$vmrs_per_chunk, "vmrs_per_chunk")
if (vmrs_per_chunk < 1L) stop("vmrs_per_chunk must be positive")
chunk_manifest <- data.frame(
    chunk_id = ceiling(tasks$task_id / vmrs_per_chunk),
    task_id = tasks$task_id,
    vmr_id = tasks$vmr_id,
    stringsAsFactors = FALSE
)

dir.create(run_dir, recursive = TRUE)
for (subdir in c("config", "logs", "results/task_rows", "results/combined")) {
    dir.create(file.path(run_dir, subdir), recursive = TRUE)
}
threshold_config <- file.path(repo_root, "config", "thresholds.yml")
invisible(file.copy(threshold_config, file.path(run_dir, "config"),
                    overwrite = FALSE))
write_tsv(tasks, file.path(run_dir, "config", "task-manifest.tsv"))
write_tsv(chunk_manifest, file.path(run_dir, "config", "chunk-manifest.tsv"))
manifest <- data.frame(
    field = c("run_id", "analysis", "run_kind", "cohort", "region",
              "started_at", "git_commit", "upstream_vmr_run_id", "vmr_set_id",
              "n_donors", "n_expected_tasks", "vmrs_per_chunk",
              "n_expected_chunks"),
    value = c(cli$run_id, "02_local_genetic_variance", "locus_geometry_scan",
              cohort, region, format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
              system2("git", c("-C", repo_root, "rev-parse", "HEAD"),
                      stdout = TRUE),
              cli$vmr_run_id, manifest_value("vmr_set_id"),
              manifest_value("n_donors"), nrow(tasks), vmrs_per_chunk,
              max(chunk_manifest$chunk_id)),
    stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(run_dir, "manifest.tsv"))
cat(normalizePath(run_dir), "\n", sep = "")
