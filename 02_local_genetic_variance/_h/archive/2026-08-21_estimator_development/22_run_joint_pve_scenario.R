#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))
source(file.path(dirname(script_path), "bslmm_pilot_functions.R"))
source(file.path(dirname(script_path), "joint_pve_functions.R"))

cli <- parse_cli(list(
    manifest = "",
    scenario_id = "",
    chunk_manifest = "",
    chunk_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = ""),
    output_root = "",
    work_root = "",
    development_bslmm_root = "",
    keep_work = "FALSE"
))
if (!nzchar(cli$manifest) || !nzchar(cli$output_root)) {
    stop("--manifest and --output_root are required")
}
manifest <- read_tsv(cli$manifest)
if (nzchar(cli$scenario_id)) {
    scenario_ids <- as_int(cli$scenario_id, "scenario_id")
} else if (nzchar(cli$chunk_manifest) && nzchar(cli$chunk_id)) {
    chunks <- read_tsv(cli$chunk_manifest)
    chunk_id <- as_int(cli$chunk_id, "chunk_id")
    scenario_ids <- as.integer(chunks$scenario_id[chunks$chunk_id == chunk_id])
    if (!length(scenario_ids)) stop("No scenarios for chunk_id=", chunk_id)
} else {
    stop("Provide --scenario_id or --chunk_manifest with --chunk_id")
}
work_root <- if (nzchar(cli$work_root)) cli$work_root else
    file.path(cli$output_root, "work")
dir.create(cli$output_root, recursive = TRUE, showWarnings = FALSE)

for (sid in scenario_ids) {
    scenario <- manifest[manifest$scenario_id == sid, , drop = FALSE]
    if (nrow(scenario) != 1L) stop("scenario_id not unique: ", sid)
    result <- run_joint_pve_features(
        scenario = scenario,
        development_bslmm_root = cli$development_bslmm_root,
        work_dir = work_root,
        keep_work = as_bool(cli$keep_work, "keep_work")
    )
    output <- file.path(cli$output_root, sprintf("scenario-%07d.tsv", sid))
    write_tsv(result, output)
    cat("Wrote", output, "\n")
}
