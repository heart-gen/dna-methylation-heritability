#!/usr/bin/env Rscript
#### 10_environmental_exploratory -- seal the run ####
##
## Sealing is not acceptance: AGENTS.md 6 makes the accepted-runs entry a human
## step, and this script deliberately does not write one.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), "10_environmental_exploratory", "_m",
                     "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

dec_f <- file.path(run_dir, "results", "environmental-decision.tsv")
if (!file.exists(dec_f)) stop("No decision file; run 04_apply_gates.R first.")
dec <- fread(dec_f)
if (startsWith(dec$decision[1], "FAIL")) {
    stop("Refusing to seal a run whose coverage gate failed: ", dec$decision[1])
}
writeLines(capture.output(sessionInfo()),
           file.path(run_dir, "results", "session-info.txt"))
append_manifest(list(dir = run_dir), list(
    decision = dec$decision[1],
    main_text_retention = dec$main_text_retention[1],
    n_eligible_exposures = as.character(dec$n_eligible_exposures[1]),
    n_tested_vmrs = as.character(dec$n_tested_vmrs[1]),
    n_axis_associations_fdr = as.character(dec$n_axis_associations_fdr[1]),
    git_commit = git_commit(), git_dirty = as.character(git_dirty()),
    sealed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
close_run(list(dir = run_dir))
message("[10] sealed ", opts$run_id, " with decision ", dec$decision[1],
        " (supplement only, always)")
