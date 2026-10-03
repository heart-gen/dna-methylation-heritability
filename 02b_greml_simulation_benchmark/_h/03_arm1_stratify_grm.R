#!/usr/bin/env Rscript
#### 02b arm 1 -- LD-stratified GRMs, one simulated sample size per task ####
##
## Port of simulation-analysis/gcta/_h/step_1b.sh + 01.stratify_LD.R:
## combine the 22 chromosome LD-score files, split SNPs into four groups at
## summary()'s Q1 / median / Q3 of ldscore_SNP, and build one GRM per group.
## Also writes a header-free copy of the phenotype file: v1 passed
## simulated.phen with its header row straight to --pheno, where GCTA read the
## header as a donor called "FID" and dropped it. Here the column order is
## checked against pheno_1..pheno_K before anything is fitted.
##
## Usage: Rscript _h/03_arm1_stratify_grm.R --run-id <id> [--task k]

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
tasks <- fread(file.path(run_dir, "config", "stratify-tasks.tsv"))
N <- tasks[task == task_id, N]
if (length(N) != 1L) stop("No stratify task ", task_id)
in_dir <- file.path(repo_root(), man[["inputs_root"]], sprintf("sim_%d_indiv", N))
w <- file.path(run_dir, "work", sprintf("sim_%d", N))
chr_files <- file.path(w, sprintf("chr%d.score.ld", seq_len(as.integer(cfg$arm1$chromosomes))))
miss <- chr_files[!file.exists(chr_files)]
if (length(miss)) stop("Missing LD score file(s): ", paste(basename(miss), collapse = ", "))
ld <- rbindlist(lapply(chr_files, fread))
if (!"ldscore_SNP" %in% names(ld)) stop("No ldscore_SNP column in GCTA LD score output")
s <- stratify_ld_quartiles(ld$ldscore_SNP)
ld[, ld_group := s$group]
fwrite(ld, file.path(w, "ldscore-combined.tsv.gz"), sep = "\t")
unlink(chr_files)

grm_list <- character(0)
for (g in seq_len(as.integer(cfg$arm1$n_ld_strata))) {
    snps <- file.path(w, sprintf("snp_group_%d.txt", g))
    writeLines(ld[ld_group == g, SNP], snps)
    pre <- file.path(w, sprintf("group_%d", g))
    r <- run_gcta(cfg$gcta$binary,
                  c("--bfile", file.path(in_dir, "plink_sim", "simulated"),
                    "--extract", snps, "--make-grm"), pre, threads = threads)
    if (r$exit != 0L || !file.exists(paste0(pre, ".grm.bin"))) {
        write_status(data.table(unit = as.character(N), status = "failed",
                                reason = paste("make-grm group", g, "exit", r$exit)),
                     file.path(run_dir, "status", sprintf("stratify_task_%d.tsv", task_id)))
        stop("make-grm failed for n=", N, " group ", g, "; see ", r$log)
    }
    grm_list <- c(grm_list, pre)
}
writeLines(grm_list, file.path(w, "multi_GRMs.txt"))

ph <- fread(file.path(in_dir, "simulated.phen"))
k <- as.integer(man[["n_phenotypes"]])
expect <- paste0("pheno_", seq_len(ncol(ph) - 2L))
if (!identical(names(ph)[-(1:2)], expect) || ncol(ph) - 2L < k) {
    stop("simulated.phen columns are not pheno_1..pheno_K in order, or fewer than ", k)
}
fwrite(ph, file.path(w, "pheno.txt"), sep = "\t", col.names = FALSE)

strata <- data.table(N = N, ld_group = seq_len(as.integer(cfg$arm1$n_ld_strata)),
                     n_snps = as.integer(table(factor(s$group, levels = 1:4))),
                     q1 = s$breaks[["q1"]], median = s$breaks[["median"]],
                     q3 = s$breaks[["q3"]])
write_atomic(strata, file.path(run_dir, "results", "stratify",
                               sprintf("strata_n%d.tsv", N)))
write_status(data.table(unit = as.character(N), status = "completed", reason = NA_character_),
             file.path(run_dir, "status", sprintf("stratify_task_%d.tsv", task_id)))
