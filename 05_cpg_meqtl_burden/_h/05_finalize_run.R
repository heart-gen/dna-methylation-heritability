#!/usr/bin/env Rscript
#### 05_cpg_meqtl_burden -- seal the run ####
##
## Sealing is not acceptance: AGENTS.md 6 makes the accepted-runs entry a human
## step, and this script deliberately does not write one.
source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), "05_cpg_meqtl_burden", "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

dec_f <- file.path(run_dir, "results", "burden-decision.tsv")
if (!file.exists(dec_f)) stop("No decision file; run 04_check_burden.R first.")
dec <- fread(dec_f)
if (startsWith(dec$decision[1], "FAIL")) {
    stop("Refusing to seal a run whose gate failed: ", dec$decision[1])
}
writeLines(capture.output(sessionInfo()),
           file.path(run_dir, "results", "session-info.txt"))
append_manifest(list(dir = run_dir), list(
    decision = dec$decision[1], git_commit = git_commit(),
    git_dirty = as.character(git_dirty()),
    sealed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
close_run(list(dir = run_dir))
message("[05] sealed ", opts$run_id, " with decision ", dec$decision[1])
