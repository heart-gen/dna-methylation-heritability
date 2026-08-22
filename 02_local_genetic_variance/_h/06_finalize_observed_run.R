#!/usr/bin/env Rscript

## Stage 06: seal a run only after the observed relative-score gate passes.
## Sealing creates provenance/checksum records and makes the run read-only; it
## does not add the run to the README's accepted-runs table.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(run_dir = "", leave_writable = "FALSE"))
if (!nzchar(cli$run_dir)) stop("--run-dir is required")
run_dir <- normalizePath(cli$run_dir)
leave_writable <- as_bool(cli$leave_writable, "leave_writable")
decision_path <- file.path(
    run_dir, "results", "combined", "observed-score-decision.tsv"
)
decision <- read_tsv(decision_path)
if (nrow(decision) != 1L ||
    !decision$decision %in% c("PASS_RELATIVE_SCORE_OBSERVED_QC",
                              "PASS_SMOKE_ONLY_NOT_ACCEPTABLE")) {
    stop("Run cannot be finalized because observed-score QC did not pass")
}

manifest_path <- file.path(run_dir, "manifest.tsv")
manifest <- read_tsv(manifest_path)
if ("finished_at" %in% manifest$field) stop("Run is already finalized")
manifest <- rbind(
    manifest,
    data.frame(
        field = c("finished_at", "observed_score_decision",
                  "absolute_pve_interpretation_allowed"),
        value = c(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                  decision$decision, "FALSE"),
        stringsAsFactors = FALSE
    )
)
write_tsv(manifest, manifest_path)
capture_session_info(file.path(run_dir, "session-info.txt"))

completed_path <- file.path(run_dir, "COMPLETED")
writeLines(c(
    paste("run_id", decision$run_id, sep = "\t"),
    paste("decision", decision$decision, sep = "\t"),
    paste("absolute_pve_authorized", "FALSE", sep = "\t"),
    paste("sealed_at", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sep = "\t")
), completed_path)

files <- list.files(run_dir, recursive = TRUE, full.names = TRUE,
                    all.files = TRUE, no.. = TRUE)
files <- files[file.info(files)$isdir %in% FALSE]
checksum_path <- file.path(run_dir, "output_checksums.tsv")
files <- setdiff(files, checksum_path)
sha256 <- vapply(files, function(path) {
    out <- system2("sha256sum", normalizePath(path), stdout = TRUE)
    if (!length(out)) stop("Could not checksum ", path)
    sub(" .*$", "", out[[1L]])
}, character(1L))
checksums <- data.frame(
    path = substring(files, nchar(run_dir) + 2L),
    sha256 = sha256,
    bytes = as.numeric(file.info(files)$size),
    stringsAsFactors = FALSE
)
checksums <- checksums[order(checksums$path), , drop = FALSE]
write_tsv(checksums, checksum_path)

if (!leave_writable) {
    status <- system2("chmod", c("-R", "a-w", run_dir))
    if (!identical(as.integer(status), 0L)) stop("Could not make run read-only")
}
cat("Finalized", decision$run_id, "with", decision$decision, "\n")
