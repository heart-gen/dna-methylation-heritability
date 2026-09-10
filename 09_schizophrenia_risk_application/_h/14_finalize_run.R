#!/usr/bin/env Rscript
#### 09_schizophrenia_risk_application -- seal the run ####
##
## Sealing is not acceptance: AGENTS.md 6 makes the accepted-runs entry a human
## step, and this script deliberately does not write one.
source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

opts <- parse_v2_args(require = "run_id")
run_dir <- file.path(repo_root(), "09_schizophrenia_risk_application", "_m",
                     "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

dec_f <- file.path(run_dir, "results", "scz-decision.tsv")
if (!file.exists(dec_f)) stop("No decision file; run 12_apply_gates.R first.")
dec <- fread(dec_f)
if (startsWith(dec$decision[1], "FAIL")) {
    stop("Refusing to seal a run whose gate failed: ", dec$decision[1])
}
writeLines(capture.output(sessionInfo()),
           file.path(run_dir, "results", "session-info.txt"))
append_manifest(list(dir = run_dir), list(
    decision = dec$decision[1],
    main_text_retention = dec$main_text_retention[1],
    n_loci_with_cpg_meqtl_support = dec$n_loci_with_cpg_meqtl_support[1],
    n_claimable_colocalizations = dec$n_claimable_colocalizations[1],
    coloc_claim_permitted = dec$coloc_claim_permitted[1],
    n_prioritized_loci = dec$n_prioritized_loci[1],
    git_commit = git_commit(), git_dirty = as.character(git_dirty()),
    sealed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
close_run(list(dir = run_dir))
message("[09] sealed ", opts$run_id, " with decision ", dec$decision[1],
        " (main-text retention: ", dec$main_text_retention[1], ")")
