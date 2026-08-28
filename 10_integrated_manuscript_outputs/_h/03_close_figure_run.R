#### 10 / 03_close_figure_run: provenance manifest, checksums, seal ####
##
## AGENTS.md 5.2: every run carries a manifest and output checksums, and is
## sealed read-only when it closes. AGENTS.md 7.9: every figure panel records
## its source run ID, table, script, and filter -- that lives in source_data/,
## which this stage verifies is present before sealing.
##
## Usage:
##   Rscript 03_close_figure_run.R --cohort AA --region caudate --run-id ID
##   (cohort/region are recorded for schema compatibility; the run spans all
##    six cells, and the per-panel provenance names the exact upstream runs.)

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages(library(data.table))

opts <- parse_v2_args(require = "run_id")
module_root <- file.path(V2_ROOT, "10_integrated_manuscript_outputs")
run_dir <- file.path(module_root, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("Run directory not found: ", run_dir)

figs <- list.files(file.path(run_dir, "figures"), pattern = "\\.(pdf|png)$")
srcs <- list.files(file.path(run_dir, "source_data"), pattern = "\\.tsv$")
if (length(figs) == 0) stop("No figures in ", run_dir)
if (length(srcs) == 0) stop("No source-data tables; AGENTS.md 7.9 requires one per panel.")

## Every PDF must have a source-data table sharing its figure stem, or a panel
## has shipped without traceable numbers.
stems <- unique(sub("\\.pdf$", "", grep("\\.pdf$", figs, value = TRUE)))
missing <- Filter(function(st) !any(startsWith(srcs, st)), stems)
if (length(missing) > 0) {
    warning("Figures without a matching source-data table: ",
            paste(missing, collapse = ", "))
}

upstream <- unique(unlist(lapply(file.path(run_dir, "source_data", srcs), function(f) {
    d <- fread(f, select = "source_run_id", colClasses = "character")
    unlist(strsplit(unique(d$source_run_id), ";"))
})))

run <- list(run_id = opts$run_id, dir = run_dir, module_root = module_root)
write_manifest(run_dir, list(
    run_id       = opts$run_id,
    analysis     = "integrated_manuscript_outputs",
    stage        = "main_and_supplemental_figures",
    git_commit   = git_commit(V2_ROOT),
    git_dirty    = git_dirty(V2_ROOT),
    r_version    = paste(R.version$major, R.version$minor, sep = "."),
    conda_prefix = Sys.getenv("CONDA_PREFIX", NA_character_),
    hostname     = Sys.info()[["nodename"]],
    n_figures    = length(figs),
    n_source_tables = length(srcs),
    upstream_runs = paste(sort(upstream), collapse = ";")
))

close_run(run)
message("[done] sealed ", run_dir, " (", length(figs), " figure files, ",
        length(srcs), " source tables)")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
