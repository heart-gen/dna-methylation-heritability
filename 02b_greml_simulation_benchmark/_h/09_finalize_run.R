#!/usr/bin/env Rscript
#### 02b -- remove regenerable intermediates and seal ####
##
## Arm 1's GRMs (up to 1.6 GB at n = 10,000) and arm 2's per-task scratch are
## regenerable from _h/, locked config and the pinned inputs, so they are
## removed before the seal rather than checksummed. What was removed is
## recorded in the manifest. The LD-score table and SNP group lists that define
## each GRM are kept.
##
## Usage: Rscript _h/09_finalize_run.R --run-id <id>

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02b_greml_simulation_benchmark", "_h")),
                 "greml_functions.R"))

opts <- parse_v2_args(require = "run_id")
run_dir <- run_dir_for(opts$run_id)
run <- list(run_id = opts$run_id, dir = run_dir)
man <- read_run_manifest(run_dir)
if (is.na(man["greml_benchmark_decision"])) stop("Gate has not run; refusing to seal")

grm <- list.files(file.path(run_dir, "work"), pattern = "\\.grm\\.(bin|N\\.bin|id)$",
                  recursive = TRUE, full.names = TRUE)
bytes <- sum(file.info(grm)$size)
unlink(grm)
scratch <- list.dirs(file.path(run_dir, "work"), recursive = TRUE)
scratch <- scratch[grepl("/(task|reml_task)_[0-9]+$", scratch)]
unlink(scratch, recursive = TRUE)
append_manifest(run, list(
    removed_before_seal = sprintf("%d GRM files (%.2f GB), %d scratch dirs; regenerable",
                                  length(grm), bytes / 1e9, length(scratch))))
close_run(run)
