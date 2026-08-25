#!/usr/bin/env Rscript
#### 03_local_snp_prediction -- seal the run ####
##
## Usage:
##   Rscript _h/05_finalize_run.R --run-id lsp-AA-caudate-20260823
##
## Records the session, checksums every output, and makes the directory
## read-only. Sealing is NOT acceptance: AGENTS.md 6 makes the accepted-runs
## entry in README.md a human step, and this script deliberately does not write
## one.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

MODULE <- "03_local_snp_prediction"
opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

comb_dir <- file.path(run_dir, "results", "combined")
dec_f <- file.path(comb_dir, "prediction-decision.tsv")
if (!file.exists(dec_f)) {
    stop("No decision file; run 04_check_prediction.R before sealing.")
}
dec <- fread(dec_f)
if (startsWith(dec$decision[1], "FAIL")) {
    stop("Refusing to seal a run whose gate failed: ", dec$decision[1])
}

writeLines(capture.output(sessionInfo()), file.path(comb_dir, "session-info.txt"))
append_manifest(list(dir = run_dir), list(
    decision = dec$decision[1],
    git_commit = git_commit(),
    git_dirty = as.character(git_dirty()),
    sealed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
))
close_run(list(dir = run_dir))
message("[03] sealed ", opts$run_id, " with decision ", dec$decision[1])
