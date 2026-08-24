#!/usr/bin/env Rscript

## Locus-geometry scan, stage B: genotype geometry for one chunk of VMRs.
##
## No phenotype is read and no model is fitted, so this costs a PLINK read and
## one crossproduct per locus. Every task writes exactly one terminal row.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
h_dir <- dirname(script_path)
source(file.path(h_dir, "00_functions.R"))
source(file.path(h_dir, "joint_pve_functions.R"))
source(file.path(h_dir, "observed_locus_io.R"))

cli <- parse_cli(list(
    run_dir = "",
    chunk_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
))
if (!nzchar(cli$run_dir) || !nzchar(cli$chunk_id)) {
    stop("--run-dir and --chunk-id are required")
}
run_dir <- normalizePath(cli$run_dir)
chunk_id <- as_int(cli$chunk_id, "chunk_id")
manifest <- read_tsv(file.path(run_dir, "manifest.tsv"))
mval <- function(field) {
    value <- manifest$value[manifest$field == field]
    if (length(value) != 1L) stop("Run manifest lacks unique field: ", field)
    as.character(value[[1L]])
}
tasks <- read_tsv(file.path(run_dir, "config", "task-manifest.tsv"))
chunks <- read_tsv(file.path(run_dir, "config", "chunk-manifest.tsv"))
wanted <- chunks$task_id[chunks$chunk_id == chunk_id]
if (!length(wanted)) stop("No tasks for chunk ", chunk_id)

repo_root <- normalizePath(file.path(h_dir, "..", ".."))
vmr_run_dir <- file.path(repo_root, "01_vmr_catalog", "_m", "runs",
                         mval("upstream_vmr_run_id"))
threshold_lines <- readLines(file.path(run_dir, "config", "thresholds.yml"),
                             warn = FALSE)
minimum_line <- grep("^[[:space:]]+min_cis_variants:", threshold_lines,
                     value = TRUE)
if (length(minimum_line) != 1L) stop("Cannot resolve min_cis_variants")
min_cis_variants <- as_int(sub(".*:[[:space:]]*", "", minimum_line),
                           "min_cis_variants")
cohort <- mval("cohort")
expected_n <- as_int(mval("n_donors"), "n_donors")

blank_row <- function(task) {
    data.frame(
        task_id = as.integer(task$task_id),
        cohort = cohort, region = mval("region"),
        chrom = as.character(task$chrom),
        start = as.integer(task$start), end = as.integer(task$end),
        n_cpgs = as.integer(task$n_cpgs),
        vmr_id = as.character(task$vmr_id),
        vmr_set_id = as.character(task$vmr_set_id),
        n = NA_integer_, num_snps = NA_integer_, snps_in_window = NA_integer_,
        p_eff = NA_real_, ld_metric = NA_real_, mean_maf = NA_real_,
        geometry_complete = FALSE, terminal_status = NA_character_,
        exclusion_reason = NA_character_, geometry_error = NA_character_,
        stringsAsFactors = FALSE
    )
}
scan_task <- function(task) {
    row <- blank_row(task)
    locus <- load_observed_locus(
        task = task, cohort = cohort, vmr_run_dir = vmr_run_dir,
        min_cis_variants = min_cis_variants, expected_n = expected_n,
        backing_tag = paste0("lgvgeo-", task$task_id)
    )
    if (!identical(locus$status, "ok")) {
        if (!is.null(locus$snps_in_window)) {
            row$snps_in_window <- as.integer(locus$snps_in_window)
        }
        row$terminal_status <- locus$status
        row$exclusion_reason <- locus$reason
        return(row)
    }
    genotype <- locus$genotype
    row$n <- nrow(genotype)
    row$num_snps <- ncol(genotype)
    row$snps_in_window <- as.integer(locus$snps_in_window)
    row$p_eff <- effective_rank_genotype(genotype)
    row$ld_metric <- adjacent_ld_metric(genotype)
    row$mean_maf <- mean(colMeans(genotype, na.rm = TRUE) / 2)
    row$geometry_complete <- all(is.finite(c(row$p_eff, row$ld_metric)))
    row$terminal_status <- if (row$geometry_complete) "completed" else
        "computational_failure"
    if (!row$geometry_complete) row$geometry_error <- "nonfinite geometry"
    row
}

for (task_id in wanted) {
    task <- tasks[tasks$task_id == task_id, , drop = FALSE]
    if (nrow(task) != 1L) stop("task_id absent or duplicated: ", task_id)
    result <- tryCatch(scan_task(task), error = function(e) {
        row <- blank_row(task)
        row$terminal_status <- "computational_failure"
        row$geometry_error <- conditionMessage(e)
        row
    })
    write_tsv(result, file.path(run_dir, "results", "task_rows",
                                sprintf("vmr-%07d.tsv", task_id)))
}
cat("Scanned", length(wanted), "loci in chunk", chunk_id, "\n")
