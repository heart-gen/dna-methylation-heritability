#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(candidate = "", repo_root = "", output = ""))
for (key in names(cli)) {
    if (!nzchar(cli[[key]])) stop("--", key, " is required")
}
candidate_path <- normalizePath(cli$candidate, mustWork = TRUE)
candidate <- read_tsv(candidate_path)
if (!"seed" %in% names(candidate)) stop("Candidate manifest lacks seed")
candidate_seed <- suppressWarnings(as.numeric(candidate$seed))
if (anyNA(candidate_seed) || anyDuplicated(candidate_seed)) {
    stop("Candidate simulation seeds must be numeric and unique")
}

roots <- file.path(normalizePath(cli$repo_root, mustWork = TRUE), c(
    "02_local_genetic_variance/_m/runs",
    "calibrated-simulation-analysis/_m/runs"
))
run_dirs <- unlist(lapply(roots[file.exists(roots)], function(root) {
    list.dirs(root, recursive = FALSE, full.names = TRUE)
}), use.names = FALSE)
config_dirs <- file.path(run_dirs, "config")
paths <- unlist(lapply(config_dirs[dir.exists(config_dirs)], function(path) {
    list.files(path, pattern = "\\.tsv$", recursive = TRUE,
               full.names = TRUE)
}), use.names = FALSE)
paths <- paths[normalizePath(paths, mustWork = TRUE) != candidate_path]

records <- list()
for (path in sort(unique(paths))) {
    header <- readLines(path, n = 1L, warn = FALSE)
    if (!length(header) || !"seed" %in% strsplit(header, "\t", fixed = TRUE)[[1L]]) {
        next
    }
    dat <- read_tsv(path)
    prior_seed <- suppressWarnings(as.numeric(dat$seed))
    prior_seed <- prior_seed[is.finite(prior_seed)]
    overlap <- intersect(candidate_seed, prior_seed)
    records[[length(records) + 1L]] <- data.frame(
        prior_manifest = normalizePath(path),
        prior_seed_rows = length(prior_seed),
        prior_unique_seeds = length(unique(prior_seed)),
        overlap_unique_seeds = length(overlap),
        stringsAsFactors = FALSE
    )
}
if (!length(records)) stop("No prior seed-bearing manifests found")
audit <- do.call(rbind, records)
audit$candidate_manifest <- candidate_path
audit$candidate_seed_rows <- length(candidate_seed)
audit$candidate_unique_seeds <- length(unique(candidate_seed))
audit$independent <- audit$overlap_unique_seeds == 0L
audit <- audit[, c("candidate_manifest", "candidate_seed_rows",
                   "candidate_unique_seeds", "prior_manifest",
                   "prior_seed_rows", "prior_unique_seeds",
                   "overlap_unique_seeds", "independent")]
write_tsv(audit, cli$output)

total_overlap <- sum(audit$overlap_unique_seeds)
cat("Audited", length(candidate_seed), "candidate seeds against",
    nrow(audit), "prior seed-bearing manifests; total overlap =",
    total_overlap, "\n")
if (total_overlap != 0L) {
    print(audit[audit$overlap_unique_seeds > 0L, , drop = FALSE])
    stop("Replacement validation seed block is not independent")
}
