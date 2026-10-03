#!/usr/bin/env Rscript
#### 02b arm 1 -- per-SNP LD scores, one (n, chromosome) per task ####
##
## Port of simulation-analysis/gcta/_h/step_1a.sh: gcta64 --ld-score-region 200
## on one chromosome of one simulated sample size.
##
## Usage: Rscript _h/02_arm1_ldscore.R --run-id <id> [--task k]

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02b_greml_simulation_benchmark", "_h")),
                 "greml_functions.R"))

opts <- parse_v2_args(require = "run_id")
run_dir <- run_dir_for(opts$run_id)
run <- list(run_id = opts$run_id, dir = run_dir)
cfg <- load_greml_config(run_dir)
man <- read_run_manifest(run_dir)
task_id <- as.integer(opts$task %||% Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
if (is.na(task_id)) stop("No --task and no SLURM_ARRAY_TASK_ID")
threads <- as.integer(Sys.getenv("V2_CPUS", "1"))
tasks <- fread(file.path(run_dir, "config", "ldscore-tasks.tsv"))
t <- tasks[task == task_id]
if (nrow(t) != 1L) stop("No ldscore task ", task_id)
bfile <- file.path(repo_root(), man[["inputs_root"]], sprintf("sim_%d_indiv", t$N),
                   "plink_sim", "simulated")
out_dir <- file.path(run_dir, "work", sprintf("sim_%d", t$N))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, sprintf("chr%d", t$chr))
r <- run_gcta(cfg$gcta$binary,
              c("--bfile", bfile, "--chr", t$chr,
                "--ld-score-region", man[["ld_score_region_kb"]]),
              out, threads = threads)
unit <- sprintf("%d:%d", t$N, t$chr)
ok <- r$exit == 0L && file.exists(paste0(out, ".score.ld"))
write_status(data.table(unit = unit, status = if (ok) "completed" else "failed",
                        reason = if (ok) NA_character_ else
                            paste("gcta_exit", r$exit, "see", basename(r$log))),
             file.path(run_dir, "status", sprintf("ldscore_task_%d.tsv", task_id)))
if (!ok) stop("LD score failed for ", unit, "; see ", r$log)
