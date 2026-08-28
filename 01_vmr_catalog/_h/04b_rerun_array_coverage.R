#### 01_vmr_catalog / 04b_rerun_qc: refresh QC tables on an accepted catalog ####
##
## Why this exists: the accepted 20260816 runs were sealed before an array probe
## universe existed on disk, so 04_turnover.R took its silent skip branch and
## every qc/array_coverage.tsv holds a placeholder note instead of the off-array
## numbers AGENTS.md 2.2 and 11 (Figure 1) require.
##
## Runs are immutable (AGENTS.md 5.2), so the fix is a NEW run that reuses the
## accepted catalog rather than an in-place edit. The VMR calls themselves are
## untouched: vmr/ is symlinked to the source run and vmr_set_id is carried
## forward unchanged, so no downstream module is invalidated.
##
## Also runs 04c_genomic_context.R, which characterizes where the corrected
## VMRs sit relative to genes -- the v2 replacement for the legacy annotation
## tree, whose outputs are keyed to the invalid legacy catalog.
##
## Usage:
##   Rscript 04b_rerun_array_coverage.R --cohort AA --region caudate \
##       --source-run-id vmrcat-AA-caudate-20260816

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages(library(data.table))

opts <- parse_v2_args(require = c("cohort", "region", "source_run_id"))
cohort <- opts$cohort; region <- opts$region
src_id <- opts$source_run_id

module_root <- file.path(V2_ROOT, "01_vmr_catalog")
src_dir <- file.path(module_root, "_m", "runs", src_id)
if (!dir.exists(src_dir)) stop("Source run not found: ", src_dir)

src_catalog <- file.path(src_dir, "vmr", "vmr_catalog.tsv")
if (!file.exists(src_catalog)) stop("Source run has no VMR catalog: ", src_catalog)
vmr_set_id <- fread(src_catalog, nrows = 1)$vmr_set_id[1]

run <- new_run("vmrcatqc", cohort, region, module_root,
               upstream = list(vmr_run_id = src_id),
               vmr_set_id = vmr_set_id,
               extra = list(
                   stage = "04_turnover_array_coverage_only",
                   rerun_reason = "array universes absent when source run was sealed"))

## Copy only the two catalog tables 04_turnover.R reads, rather than symlinking
## vmr/. A symlinked directory would put the sealed source run inside the reach
## of close_run()'s recursive chmod, and copying makes it explicit that this run
## re-derives nothing: the catalog arrives verbatim, checksums included.
dir.create(file.path(run$dir, "vmr"), showWarnings = FALSE)
for (f in c("vmr_catalog.tsv", "cpg_vmr_membership.tsv")) {
    ok <- file.copy(file.path(src_dir, "vmr", f), file.path(run$dir, "vmr", f))
    if (!ok) stop("Failed to copy ", f, " from ", src_id)
    Sys.chmod(file.path(run$dir, "vmr", f), mode = "0644")
}

## prepare_summary.tsv files feed the technical-QC table; copy the tree shape.
prep <- list.files(src_dir, pattern = "^prepare_summary\\.tsv$",
                   recursive = TRUE, full.names = TRUE)
for (f in prep) {
    rel <- sub(paste0("^", src_dir, "/"), "", f)
    dest <- file.path(run$dir, rel)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(f, dest, overwrite = TRUE)
}
message("[stage] copied catalog tables and ", length(prep), " prepare summaries")

run_stage <- function(script) {
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c(shQuote(file.path(module_root, "_h", script)),
                        "--cohort", cohort, "--region", region,
                        "--run-id", run$run_id))
    if (status != 0L) stop(script, " failed with exit status ", status)
}
run_stage("04_turnover.R")
## Genomic context is descriptive characterization of the same catalog, so it
## belongs in the same refreshed QC run rather than a third run directory.
run_stage("04c_genomic_context.R")

cov <- fread(file.path(run$dir, "qc", "array_coverage.tsv"))
if (!"array_platform" %in% names(cov) || anyNA(cov$frac_vmr_cpgs_off_array)) {
    stop("array_coverage.tsv is still missing off-array numbers; refusing to seal.")
}
print(cov[, .(array_platform, frac_vmr_cpgs_off_array, frac_vmrs_invisible_to_array)])

close_run(run)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
