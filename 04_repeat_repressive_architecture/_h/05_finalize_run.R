#!/usr/bin/env Rscript
#### 04_repeat_repressive_architecture -- seal the cells ####
##
## Sealing is not acceptance: AGENTS.md 6 makes the accepted-runs entry a human
## step, and this script deliberately does not write one.
##
## Unlike modules 03 and 05, this module's gate is CROSS-REGION: 03_apply_gates.R
## writes one claims table into the FIRST run's results directory, and it is that
## table -- not any per-cell decision -- that says what the module is allowed to
## conclude. So every cell is sealed together, and none is sealed unless the
## claims table exists.
##
##   Rscript _h/05_finalize_run.R --cohort AA --run-ids id1,id2,id3

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
suppressPackageStartupMessages(library(data.table))

MODULE <- "04_repeat_repressive_architecture"
opts <- parse_v2_args(require = c("cohort", "run_ids"))
run_ids <- trimws(strsplit(opts$run_ids, ",")[[1]])

run_dirs <- file.path(repo_root(), MODULE, "_m", "runs", run_ids)
missing <- run_ids[!dir.exists(run_dirs)]
if (length(missing) > 0) stop("No such run(s): ", paste(missing, collapse = ", "))

claims_f <- file.path(run_dirs[1], "results", "interpretation-claims.tsv")
if (!file.exists(claims_f)) {
    stop("No interpretation-claims.tsv under ", run_ids[1],
         "; run 03_apply_gates.R across all regions before finalizing.")
}
claims <- fread(claims_f)

## The decision records what the gates permit, not whether the code ran. A run
## in which no outcome survives is a legitimate, sealable negative result -- the
## failure mode this module must avoid is sealing without the gates having been
## applied at all, which the check above catches.
n_supported <- sum(!startsWith(claims$permitted_claim, "not supported"))
decision <- sprintf("GATES_APPLIED_%d_OF_%d_OUTCOMES_SUPPORTED",
                    n_supported, nrow(claims))

for (rd in run_dirs) {
    smoke <- {
        m <- fread(file.path(rd, "manifest.tsv"))
        identical(m$value[m$field == "smoke_run"][1], "TRUE")
    }
    writeLines(capture.output(sessionInfo()),
               file.path(rd, "results", "session-info.txt"))
    append_manifest(list(dir = rd), list(
        decision = if (smoke) paste0(decision, "_SMOKE_ONLY_NOT_ACCEPTABLE")
                   else decision,
        gate_run_id = run_ids[1],
        git_commit = git_commit(),
        git_dirty = as.character(git_dirty()),
        sealed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
    close_run(list(dir = rd))
}
message("[04] sealed ", length(run_dirs), " cell(s) with decision ", decision)
