#!/usr/bin/env Rscript
#### 06_partitioned_heritability -- seal the run ####
##
## Sealing is not acceptance: AGENTS.md 6 makes the accepted-runs entry a human
## step, and this script deliberately does not write one.
source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), "06_partitioned_heritability", "_m", "runs",
                     opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

dec_f <- file.path(run_dir, "results", "partitioned-h2-decision.tsv")
if (!file.exists(dec_f)) stop("No decision file; run 07_fdr_and_gates.R first.")
dec <- fread(dec_f)
if (startsWith(dec$decision[1], "FAIL")) {
    stop("Refusing to seal a run whose gate failed: ", dec$decision[1])
}

writeLines(capture.output(sessionInfo()),
           file.path(run_dir, "results", "session-info.txt"))
append_manifest(list(dir = run_dir), list(
    decision = dec$decision[1],
    n_brain_significant = dec$n_brain_significant[1],
    n_control_significant = dec$n_control_significant[1],
    sldsc_supports_brain_enrichment = dec$sldsc_supports_brain_enrichment[1],
    git_commit = git_commit(),
    git_dirty = as.character(git_dirty()),
    sealed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
close_run(list(dir = run_dir))
message("[06] sealed ", opts$run_id, " with decision ", dec$decision[1])
