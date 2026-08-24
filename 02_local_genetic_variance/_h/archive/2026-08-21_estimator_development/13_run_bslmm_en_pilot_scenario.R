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
    calibration_model = "",
    expected_calibration_sha256 = "",
    output_root = "",
    work_root = "",
    keep_work = "FALSE"
))
if (!nzchar(cli$manifest) || !nzchar(cli$output_root)) {
    stop("--manifest and --output_root are required")
}
if (!nzchar(cli$calibration_model) || !nzchar(cli$expected_calibration_sha256)) {
    stop("--calibration_model and --expected_calibration_sha256 are required")
}
verify_calibration_sha256(cli$calibration_model, cli$expected_calibration_sha256)
model <- readRDS(cli$calibration_model)
manifest <- read_tsv(cli$manifest)

scenario_ids <- integer()
if (nzchar(cli$scenario_id)) {
    scenario_ids <- as_int(cli$scenario_id, "scenario_id")
} else if (nzchar(cli$chunk_manifest) && nzchar(cli$chunk_id)) {
    chunk_id <- as_int(cli$chunk_id, "chunk_id")
    chunks <- read_tsv(cli$chunk_manifest)
    if (!all(c("chunk_id", "scenario_id") %in% names(chunks))) {
        stop("chunk manifest must contain chunk_id and scenario_id")
    }
    scenario_ids <- as.integer(chunks$scenario_id[chunks$chunk_id == chunk_id])
    if (!length(scenario_ids)) stop("No scenarios for chunk_id=", chunk_id)
} else {
    stop("Provide --scenario_id or --chunk_manifest with --chunk_id / SLURM_ARRAY_TASK_ID")
}

work_root <- if (nzchar(cli$work_root)) {
    cli$work_root
} else {
    file.path(cli$output_root, "work")
}
keep_work <- as_bool(cli$keep_work, "keep_work")
dir.create(cli$output_root, recursive = TRUE, showWarnings = FALSE)

for (sid in scenario_ids) {
    scenario <- manifest[manifest$scenario_id == sid, , drop = FALSE]
    if (nrow(scenario) != 1L) stop("scenario_id not unique in manifest: ", sid)
    result <- run_paired_en_bslmm(
        scenario = scenario,
        calibration_model = model,
        gemma_bin = scenario$gemma_bin[[1L]],
        bslmm_mode = as.integer(scenario$bslmm_mode[[1L]]),
        burn_in = as.integer(scenario$bslmm_burn_in[[1L]]),
        sampling = as.integer(scenario$bslmm_sampling[[1L]]),
        rpace = as.integer(scenario$bslmm_rpace[[1L]]),
        work_dir = work_root,
        keep_work = keep_work
    )
    output <- file.path(
        cli$output_root, sprintf("scenario-%07d.tsv", sid)
    )
    write_tsv(result, output)
    cat("Wrote", output, "\n")
}
