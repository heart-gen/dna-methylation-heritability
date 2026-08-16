#### 01_vmr_catalog / 05_close_run: checksum outputs and seal the run ####
##
## The last thing a run does. Records the `vmr_set_id` and the observed sample
## count in the manifest, checksums every output, and makes the directory
## read-only so it cannot be edited in place (AGENTS.md §5.2, §9).
##
## Usage:
##   Rscript 05_close_run.R --cohort AA --region caudate --run-id ID

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages(library(data.table))

opts <- parse_v2_args(require = c("cohort", "region", "run_id"))

module_root <- file.path(V2_ROOT, "01_vmr_catalog")
run_dir <- file.path(module_root, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("Run directory not found: ", run_dir)

summary_file <- file.path(run_dir, "vmr", "summarize_summary.tsv")
if (!file.exists(summary_file)) {
    stop("Cannot close a run that has not been summarized: ", summary_file,
         "\n  Run 02_summarize.R first.")
}
s <- fread(summary_file)
field <- function(f) {
    i <- match(f, s$field)
    if (is.na(i)) NA_character_ else s$value[i]
}

run <- list(run_id = opts$run_id, dir = run_dir, module_root = module_root)

append_manifest(run, list(
    vmr_set_id     = field("vmr_set_id"),
    n_vmrs         = field("n_vmrs"),
    n_donors       = field("n_donors"),
    n_cpgs_in_vmrs = field("n_cpgs_in_vmrs"),
    donor_checksum = field("donor_checksum"),
    bsseq_version  = field("bsseq_version"),
    chromosomes    = field("chromosomes")
))

close_run(run)

message("[close] run ", opts$run_id, " is sealed.\n",
        "  vmr_set_id: ", field("vmr_set_id"), "\n",
        "  Record it in 01_vmr_catalog/README.md with its acceptance gate ",
        "before any downstream module consumes it (AGENTS.md 6).")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
