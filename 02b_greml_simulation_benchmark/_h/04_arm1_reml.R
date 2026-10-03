#!/usr/bin/env Rscript
#### 02b arm 1 -- GREML-LDMS on a chunk of simulated phenotypes ####
##
## Port of simulation-analysis/gcta/_h/step_1c.sh (gcta64 --reml --mgrm,
## --mpheno k, no covariates), run once per REML mode in config. The v1
## setting is the constrained mode; the unconstrained mode is added so arm 1
## and arm 2 are read on the same estimator.
##
## Usage: Rscript _h/04_arm1_reml.R --run-id <id> [--task k]

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
tasks <- fread(file.path(run_dir, "config", "reml-tasks.tsv"))
t <- tasks[task == task_id]
if (nrow(t) != 1L) stop("No reml task ", task_id)
modes <- fread(file.path(run_dir, "config", "reml-modes.tsv"))
w <- file.path(run_dir, "work", sprintf("sim_%d", t$N))
mgrm <- file.path(w, "multi_GRMs.txt")
if (!file.exists(mgrm)) stop("No GRM list for n=", t$N, "; did stratify finish?")
tmp <- file.path(w, sprintf("reml_task_%d", task_id))
dir.create(tmp, showWarnings = FALSE)
rows <- list(); st <- list()
for (k in seq(t$pheno_from, t$pheno_to)) for (i in seq_len(nrow(modes))) {
    m <- modes[i]
    unit <- sprintf("%d:pheno_%d:%s", t$N, k, m$reml_mode)
    out <- file.path(tmp, sprintf("pheno_%d_%s", k, m$reml_mode))
    margs <- if (nzchar(m$args)) strsplit(m$args, " ", fixed = TRUE)[[1]] else character(0)
    r <- run_gcta(cfg$gcta$binary,
                  c("--reml", "--mgrm", mgrm, "--pheno", file.path(w, "pheno.txt"),
                    "--mpheno", k, "--reml-maxit", cfg$gcta$reml_maxit, margs),
                  out, threads = threads)
    cl <- classify_reml(r$exit, paste0(out, ".hsq"), r$log,
                        unlist(cfg$gcta$estimation_failure_patterns))
    st[[length(st) + 1L]] <- data.table(unit = unit, status = cl$status, reason = cl$reason)
    if (cl$status == "completed") {
        rows[[length(rows) + 1L]] <- cbind(
            data.table(arm = man[["arm"]], N = t$N, phenotype_id = paste0("pheno_", k),
                       reml_mode = m$reml_mode, constrained = m$constrained),
            parse_hsq(paste0(out, ".hsq")))
    }
    unlink(list.files(tmp, pattern = paste0("^", basename(out), "\\."), full.names = TRUE))
}
res <- rbindlist(rows, fill = TRUE)
write_atomic(res, file.path(run_dir, "results", "reml", sprintf("task_%d.tsv", task_id)))
write_status(rbindlist(st), file.path(run_dir, "status", sprintf("reml_task_%d.tsv", task_id)))
unlink(tmp, recursive = TRUE)
failed <- rbindlist(st)[status == "failed"]
if (nrow(failed)) stop(nrow(failed), " REML call(s) failed computationally; first: ",
                       failed$reason[1])
