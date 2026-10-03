#!/usr/bin/env Rscript
#### 02b -- collate accepted runs into _m/combined/ ####
##
## Reads only runs listed under "Accepted runs" in the module README. With
## --allow-unaccepted it collates the newest sealed run of each cell instead
## and writes -UNACCEPTED filenames, which are gitignored and not citable.
##
## Usage: Rscript _h/10_collate.R [--allow-unaccepted TRUE]

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02b_greml_simulation_benchmark", "_h")),
                 "greml_functions.R"))

opts <- parse_v2_args(require = character(0))
allow <- identical(toupper(opts$allow_unaccepted %||% "FALSE"), "TRUE")
module <- "02b_greml_simulation_benchmark"
cfg <- load_config("greml_benchmark")
cells <- rbind(data.table(cohort = cfg$arm1$run_cohort, region = cfg$arm1$run_region),
               data.table(cohort = cfg$arm2$cohort, region = unlist(cfg$arm2$regions)))
ids <- vapply(seq_len(nrow(cells)), function(i) {
    r <- tryCatch(require_accepted_upstream(module, cells$cohort[i], cells$region[i]),
                  error = function(e) NULL)
    if (!is.null(r)) return(r$run_id)
    if (!allow) stop("No accepted 02b run for ", cells$cohort[i], " x ", cells$region[i])
    runs <- list.files(file.path(repo_root(), module, "_m", "runs"),
                       pattern = sprintf("^greml-%s-%s-", cells$cohort[i], cells$region[i]))
    sealed <- runs[vapply(runs, function(x) run_is_sealed(
        fread(file.path(run_dir_for(x), "manifest.tsv"), colClasses = "character")),
        logical(1))]
    if (!length(sealed)) return(NA_character_)
    sort(sealed, decreasing = TRUE)[1]
}, character(1))
cells[, run_id := ids]
cells <- cells[!is.na(run_id)]
suffix <- if (allow) "-UNACCEPTED" else ""
out <- file.path(repo_root(), module, "_m", "combined")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
grab <- function(f) rbindlist(lapply(cells$run_id, function(id) {
    d <- fread(file.path(run_dir_for(id), "summary", f))
    d[, source_run_id := id]
}), fill = TRUE)
write_atomic(grab("recovery-metrics.tsv"),
             file.path(out, paste0("greml-benchmark-recovery-metrics", suffix, ".tsv")))
write_atomic(grab("spearman.tsv"),
             file.path(out, paste0("greml-benchmark-spearman", suffix, ".tsv")))
dec <- rbindlist(lapply(cells$run_id, function(id) {
    m <- read_run_manifest(run_dir_for(id))
    data.table(run_id = id, arm = m[["arm"]], cohort = m[["cohort"]], region = m[["region"]],
               decision = m[["greml_benchmark_decision"]],
               gate_criteria_passed = m[["gate_criteria_passed"]],
               gcta_version = m[["gcta_version"]], smoke_run = m[["smoke_run"]],
               citable = !allow)
}))
write_atomic(dec, file.path(out, paste0("greml-benchmark-decisions", suffix, ".tsv")))
print(dec)
