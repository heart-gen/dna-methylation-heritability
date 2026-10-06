#!/usr/bin/env Rscript
#### 02b arm 1 -- pin the simulated inputs ####
##
## Checksums every simulated input the arm reads (AGENTS.md 9). Separate from
## 00_new_run.R because the PLINK files at n = 5,000 and 10,000 are large enough
## that hashing them belongs on a compute node, not on the submit host.
##
## Usage: Rscript _h/01_arm1_inputs.R --run-id <id>

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
source(file.path(Sys.getenv("V2_RUN_CODE",
                            file.path(repo_root(), "02b_greml_simulation_benchmark", "_h")),
                 "greml_functions.R"))

opts <- parse_v2_args(require = "run_id")
run_dir <- run_dir_for(opts$run_id)
run <- list(run_id = opts$run_id, dir = run_dir)
cfg <- load_greml_config(run_dir)
man <- read_run_manifest(run_dir)
sizes <- as.integer(strsplit(man[["sample_sizes"]], ",")[[1]])
root_in <- file.path(repo_root(), man[["inputs_root"]])
files <- unlist(lapply(sizes, function(N) {
    d <- file.path(root_in, sprintf("sim_%d_indiv", N))
    c(file.path(d, "plink_sim", paste0("simulated.", c("bed", "bim", "fam"))),
      file.path(d, "simulated.phen"), file.path(d, "snp_phenotype_mapping.tsv"))
}))
sums <- data.table(
    file = sub(paste0("^", repo_root(), "/"), "", files),
    resolved = normalizePath(files),
    bytes = file.info(files)$size,
    sha256 = vapply(files, file_sha256, character(1)))
write_atomic(sums, file.path(run_dir, "inputs", "input_checksums.tsv"))
append_manifest(run, list(n_input_files = nrow(sums),
                          input_checksums = "inputs/input_checksums.tsv"))
message("[02b] pinned ", nrow(sums), " arm 1 input files")
