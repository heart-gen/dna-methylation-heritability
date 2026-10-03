#!/usr/bin/env Rscript
#### 02c -- collate accepted runs into _m/combined/ ####
##
## Reads only runs listed under "Accepted runs" in the module README. With
## --allow-unaccepted TRUE it takes the newest sealed run per region and writes
## -UNACCEPTED filenames (gitignored, not citable). Region results are stacked,
## never compared at level across regions (AGENTS.md 7.6).
##
## Usage: Rscript _h/06_collate.R [--allow-unaccepted TRUE]

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02c_cis_greml_sensitivity", "_h")),
                 "cgs_functions.R"))
opts <- parse_v2_args(require = character(0))
allow <- identical(toupper(opts$allow_unaccepted %||% "FALSE"), "TRUE")
cfg <- load_config("cis_greml_sensitivity")
ids <- vapply(unlist(cfg$regions), function(r) {
    a <- tryCatch(require_accepted_upstream(MODULE, cfg$cohort, r), error = function(e) NULL)
    if (!is.null(a)) return(a$run_id)
    if (!allow) stop("No accepted 02c run for ", cfg$cohort, " x ", r)
    runs <- list.files(file.path(repo_root(), MODULE, "_m", "runs"),
                       pattern = sprintf("^cgs-%s-%s-[0-9]{8}(-[a-z])?$", cfg$cohort, r))
    sealed <- runs[vapply(runs, function(x) run_is_sealed(fread(
        file.path(run_dir_for(x), "manifest.tsv"), colClasses = "character")), logical(1))]
    if (length(sealed)) sort(sealed, decreasing = TRUE)[1] else NA_character_
}, character(1))
ids <- ids[!is.na(ids)]
sfx <- if (allow) "-UNACCEPTED" else ""
out <- file.path(repo_root(), MODULE, "_m", "combined")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
for (f in c("existence", "ordering", "decile-profile", "convergence-by-decile")) {
    d <- rbindlist(lapply(ids, function(id) fread(file.path(run_dir_for(id), "summary",
                                                            paste0(f, ".tsv")))[, source_run_id := id]),
                   fill = TRUE)
    d[, citable := !allow]
    write_atomic(d, file.path(out, paste0("cis-greml-", f, "-", cfg$cohort, sfx, ".tsv")))
}
print(ids)
