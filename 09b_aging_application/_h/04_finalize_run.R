#!/usr/bin/env Rscript
#### 09b_aging_application -- seal the run ####
##
## Sealing is not acceptance: AGENTS.md 6 makes the accepted-runs entry a human
## step, and this script deliberately does not write one. The donor x VMR
## checkpoint and the bootstrap draws stay in the run: stage 05 bootstraps
## against the checkpoint, and the draws are half of the headline's variance.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), "09b_aging_application", "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

dec_f <- file.path(run_dir, "results", "aging-decision.tsv")
if (!file.exists(dec_f)) stop("No decision file; run 03_apply_gates.R first.")
dec <- fread(dec_f)
if (startsWith(dec$decision[1], "FAIL")) {
    stop("Refusing to seal a run whose coverage gate failed: ", dec$decision[1])
}
writeLines(capture.output(sessionInfo()),
           file.path(run_dir, "results", "session-info.txt"))
append_manifest(list(dir = run_dir), list(
    decision = dec$decision[1],
    region_reading = dec$region_reading[1],
    primary_p = format(dec$primary_p[1], digits = 4),
    git_commit = git_commit(), git_dirty = as.character(git_dirty()),
    sealed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
close_run(list(dir = run_dir))
message("[09b] sealed ", opts$run_id, ": ", dec$decision[1], ", ",
        dec$region_reading[1])
