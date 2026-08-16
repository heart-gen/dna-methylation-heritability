#### 01_vmr_catalog / 00_new_run: create the immutable run directory ####
##
## Called once per cohort x region before any chromosome fans out, so every
## array task writes into one run whose provenance was recorded up front.
## Prints the run ID on stdout for the launcher to capture.
##
## Usage:
##   Rscript 00_new_run.R --cohort AA --region caudate [--allow-unlocked]

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

opts <- parse_v2_args(require = c("cohort", "region"))

assert_locked(
    list(cohorts = load_config("cohorts"),
         thresholds = load_config("thresholds"),
         covariates_vmr = list(pi_locked = load_config("covariates")$vmr_calling$pi_locked)),
    allow_unlocked = opts$allow_unlocked)

run <- new_run(
    module = "vmrcat",
    cohort = opts$cohort,
    region = opts$region,
    module_root = file.path(V2_ROOT, "01_vmr_catalog"),
    extra = list(
        smoke_run = opts$allow_unlocked,
        chromosome_policy = paste(chrom_order(), collapse = ","),
        cis_window_bp = load_config("thresholds")$cis$window_bp
    ))

## The launcher reads this line.
cat(run$run_id, "\n", sep = "")
