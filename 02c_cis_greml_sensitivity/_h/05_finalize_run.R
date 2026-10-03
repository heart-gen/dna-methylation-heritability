#!/usr/bin/env Rscript
#### 02c -- remove scratch and seal ####
## Usage: Rscript _h/05_finalize_run.R --run-id <id>

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02c_cis_greml_sensitivity", "_h")),
                 "cgs_functions.R"))

opts <- parse_v2_args(require = "run_id")
run_dir <- run_dir_for(opts$run_id)
run <- list(run_id = opts$run_id, dir = run_dir)
cfg <- load_cgs_config(run_dir)
man <- read_run_manifest(run_dir)
smoke <- identical(man[["smoke_run"]], "TRUE")
if (is.na(man["cis_greml_sensitivity_decision"])) stop("Gate has not run; refusing to seal")
unlink(file.path(run_dir, "work"), recursive = TRUE)
close_run(run)
