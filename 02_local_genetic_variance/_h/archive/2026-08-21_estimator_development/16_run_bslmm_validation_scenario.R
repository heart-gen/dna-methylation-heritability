#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))
source(file.path(dirname(script_path), "bslmm_pilot_functions.R"))

cli <- parse_cli(list(
    manifest = "",
    scenario_id = "",
    chunk_manifest = "",
    chunk_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = ""),
    output_root = "",
    work_root = "",
    keep_work = "FALSE"
))
if (!nzchar(cli$manifest) || !nzchar(cli$output_root)) {
    stop("--manifest and --output_root are required")
}
manifest <- read_tsv(cli$manifest)
scenario_ids <- integer()
if (nzchar(cli$scenario_id)) {
    scenario_ids <- as_int(cli$scenario_id, "scenario_id")
} else if (nzchar(cli$chunk_manifest) && nzchar(cli$chunk_id)) {
    chunk_id <- as_int(cli$chunk_id, "chunk_id")
    chunks <- read_tsv(cli$chunk_manifest)
    scenario_ids <- as.integer(chunks$scenario_id[chunks$chunk_id == chunk_id])
    if (!length(scenario_ids)) stop("No scenarios for chunk_id=", chunk_id)
} else {
    stop("Provide --scenario_id or --chunk_manifest with chunk id")
}

work_root <- if (nzchar(cli$work_root)) cli$work_root else file.path(cli$output_root, "work")
keep_work <- as_bool(cli$keep_work, "keep_work")
dir.create(cli$output_root, recursive = TRUE, showWarnings = FALSE)

for (sid in scenario_ids) {
    scenario <- manifest[manifest$scenario_id == sid, , drop = FALSE]
    if (nrow(scenario) != 1L) stop("scenario_id not unique: ", sid)
    result <- run_bslmm_only(
        scenario = scenario,
        gemma_bin = scenario$gemma_bin[[1L]],
        bslmm_mode = as.integer(scenario$bslmm_mode[[1L]]),
        burn_in = as.integer(scenario$bslmm_burn_in[[1L]]),
        sampling = as.integer(scenario$bslmm_sampling[[1L]]),
        rpace = as.integer(scenario$bslmm_rpace[[1L]]),
        work_dir = work_root,
        keep_work = keep_work,
        positive_signal_rule = scenario$positive_signal_rule[[1L]]
    )
    write_tsv(result, file.path(cli$output_root, sprintf("scenario-%07d.tsv", sid)))
    cat("Wrote scenario", sid, "\n")
}
