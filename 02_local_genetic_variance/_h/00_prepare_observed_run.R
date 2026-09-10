#!/usr/bin/env Rscript

## Stage 00: open one immutable observed local-genetic-control run.
## This stage performs no estimation. It resolves and freezes the accepted VMR
## catalog, final joint model, task universe, settings, and checksums.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
h_dir <- dirname(script_path)
source(file.path(h_dir, "00_functions.R"))

repo_root <- normalizePath(file.path(h_dir, "..", ".."))
module_root <- file.path(repo_root, "02_local_genetic_variance")
cli <- parse_cli(list(
    run_id = "",
    cohort = "",
    region = "",
    vmr_run_id = "",
    model_run_id = "lgv-joint-pve-train-20260820",
    model_sha256 =
        "9f26c3273746fda85d9bbf21e224857db9a1ad79a521582a12f241854c03223a",
    smoke_n = "0",
    vmrs_per_chunk = "5",
    runs_root = ""
))

required_cli <- c("run_id", "cohort", "region", "vmr_run_id")
missing_cli <- required_cli[!vapply(cli[required_cli], nzchar, logical(1L))]
if (length(missing_cli)) {
    stop("Required argument(s): ", paste(missing_cli, collapse = ", "))
}
cohort <- cli$cohort
region <- tolower(cli$region)

## `cohort` is the CELL token. It is either a discovery arm (AA,
## all_individuals) or an estimation cell (all_individuals.AA,
## all_individuals.EA) -- a donor group modelled on a locus set discovered
## somewhere else. config/cohorts.yml is the single source of that list; this
## stage used to carry a two-element literal, which is what made discovery and
## estimation the same population by construction.
## 00_shared/load.R defines repo_root() as a FUNCTION, and it sources its
## members into the global environment. This script bound `repo_root` to a PATH
## twenty lines above and uses it in every file.path() below, so sourcing here
## silently replaces the path with the closure. Save it across the source; the
## two meanings of the name are the hazard, not the source itself.
.repo_root_path <- repo_root
source(file.path(.repo_root_path, "00_shared", "load.R"))
repo_root <- .repo_root_path
rm(.repo_root_path)

validate_cohort_region(cohort, region, root = repo_root)
cell <- parse_cell(cohort, root = repo_root)
catalog_cohort <- cell$catalog_cohort
estimation_group <- cell$estimation_group

## The cell token may carry a dot, so it is quoted into the pattern rather than
## interpolated raw -- otherwise "all_individuals.AA" would also match
## "all_individualsXAA".
expected_run_pattern <- paste0(
    "^lgv-", gsub(".", "\\.", cohort, fixed = TRUE), "-", region,
    "-[0-9]{8}[a-z]?$"
)
if (!grepl(expected_run_pattern, cli$run_id)) {
    stop("run_id must match ", expected_run_pattern)
}

runs_root <- if (nzchar(cli$runs_root)) cli$runs_root else
    file.path(module_root, "_m", "runs")
dir.create(runs_root, recursive = TRUE, showWarnings = FALSE)
run_dir <- file.path(runs_root, cli$run_id)
if (file.exists(run_dir)) stop("Run directory already exists: ", run_dir)

## An estimation cell is materialized by 01b_estimation_cells, which writes a
## run satisfying the same directory contract as a Module 01 run (vmr/,
## covs/, plink_format/). Route by the run-ID prefix rather than by the cohort,
## so a cell can be smoke-tested against either module while the two coexist.
upstream_module <- if (startsWith(cli$vmr_run_id, "estcell-")) {
    "01b_estimation_cells"
} else {
    "01_vmr_catalog"
}
vmr_run_dir <- file.path(
    repo_root, upstream_module, "_m", "runs", cli$vmr_run_id
)
if (cell$is_estimation_cell && !identical(upstream_module,
                                          "01b_estimation_cells")) {
    stop("Cell '", cohort, "' must consume an 01b_estimation_cells run. ",
         "A Module 01 run holds the POOLED donors, so Module 02 would silently ",
         "estimate in the wrong donor set. Got: ", cli$vmr_run_id)
}
if (!cell$is_estimation_cell && identical(upstream_module,
                                          "01b_estimation_cells")) {
    stop("Arm '", cohort, "' must consume a 01_vmr_catalog run, not ",
         cli$vmr_run_id)
}
vmr_manifest_path <- file.path(vmr_run_dir, "manifest.tsv")
vmr_catalog_path <- file.path(vmr_run_dir, "vmr", "vmr_catalog.tsv")
for (path in c(vmr_manifest_path, vmr_catalog_path)) {
    if (!file.exists(path)) stop("Required Module 01 input is missing: ", path)
}

## Module 01 predates the shared `Accepted runs` parser. Its README is still
## the record of acceptance, so require the exact run row and its locked gate.
## 01b_estimation_cells uses the same README convention and the same gate
## wording, so one check covers both upstreams.
vmr_readme <- readLines(file.path(repo_root, upstream_module, "README.md"),
                        warn = FALSE)
acceptance_row <- vmr_readme[grepl(
    paste0("| `", cli$vmr_run_id, "` |"), vmr_readme, fixed = TRUE
)]
if (length(acceptance_row) != 1L ||
    !grepl("all five pass", acceptance_row, fixed = TRUE)) {
    stop("Upstream ", upstream_module, " run is not recorded as passing all ",
         "five gates: ", cli$vmr_run_id)
}

vmr_manifest <- read_tsv(vmr_manifest_path)
manifest_value <- function(field) {
    value <- vmr_manifest$value[vmr_manifest$field == field]
    if (length(value) != 1L) stop("Module 01 manifest lacks unique field: ", field)
    as.character(value[[1L]])
}
## The upstream run's own cohort field is the CELL it carries donors for: a
## Module 01 run is its arm, an 01b run is its cell. Either way it must equal
## the cell being estimated -- this is what pins Module 02 to the right donors.
if (!identical(manifest_value("cohort"), cohort) ||
    !identical(tolower(manifest_value("region")), region)) {
    stop("Requested cohort/region does not match the ", upstream_module,
         " manifest")
}
## And, separately, the loci must come from the cell's DISCOVERY arm. For an
## arm these are the same statement; for an estimation cell they are not, and
## keeping them separate is what lets AA and EA share one pooled locus set
## (AGENTS.md 7.7).
if (cell$is_estimation_cell) {
    if (!identical(manifest_value("catalog_cohort"), catalog_cohort)) {
        stop("Cell '", cohort, "' discovers on '", catalog_cohort,
             "' but its 01b run was built on '",
             manifest_value("catalog_cohort"), "'")
    }
    if (!identical(manifest_value("estimation_group"), estimation_group)) {
        stop("Cell '", cohort, "' estimates in group '", estimation_group,
             "' but its 01b run carries '",
             manifest_value("estimation_group"), "'")
    }
}

## The covariate files inside the run carry the DISCOVERY arm's prefix, because
## that arm's Module 01 run wrote them. locus_io.R cannot infer this from the
## cell token, so it is resolved once here and travels in the manifest.
covar_prefix <- if (cell$is_estimation_cell) {
    manifest_value("covar_prefix")
} else if (identical(cohort, "AA")) {
    "TOPMed_LIBD.AA"
} else {
    "TOPMed_LIBD"
}
if (!identical(toupper(manifest_value("smoke_run")), "FALSE")) {
    stop("Observed production must start from a non-smoke Module 01 run")
}

model_dir <- file.path(
    module_root, "_m", "runs", cli$model_run_id, "combined"
)
model_path <- file.path(model_dir, "joint-pve-calibrator.rds")
development_features <- file.path(model_dir, "development-features.tsv")
for (path in c(model_path, development_features)) {
    if (!file.exists(path)) stop("Frozen joint-model input is missing: ", path)
}
sha256_file <- function(path) {
    out <- system2("sha256sum", normalizePath(path), stdout = TRUE)
    if (!length(out)) stop("Could not checksum ", path)
    tolower(sub(" .*$", "", out[[1L]]))
}
observed_model_sha <- sha256_file(model_path)
if (!identical(observed_model_sha, tolower(cli$model_sha256))) {
    stop("Frozen joint-model checksum mismatch\n  expected ", cli$model_sha256,
         "\n  observed ", observed_model_sha)
}

tasks <- read_tsv(vmr_catalog_path)
required_task <- c("chr", "start", "end", "n", "vmr_id", "vmr_set_id")
missing_task <- setdiff(required_task, names(tasks))
if (length(missing_task)) {
    stop("VMR catalog lacks: ", paste(missing_task, collapse = ", "))
}
if (anyDuplicated(tasks$vmr_id)) stop("Duplicate vmr_id in VMR catalog")
if (length(unique(tasks$vmr_set_id)) != 1L ||
    !identical(unique(tasks$vmr_set_id), manifest_value("vmr_set_id"))) {
    stop("VMR catalog and manifest vmr_set_id differ")
}
tasks <- data.frame(
    task_id = seq_len(nrow(tasks)),
    cohort = cohort,
    region = region,
    chrom = as.character(tasks$chr),
    start = as.integer(tasks$start),
    end = as.integer(tasks$end),
    n_cpgs = as.integer(tasks$n),
    vmr_id = as.character(tasks$vmr_id),
    vmr_set_id = as.character(tasks$vmr_set_id),
    stringsAsFactors = FALSE
)
smoke_n <- as_int(cli$smoke_n, "smoke_n")
if (smoke_n < 0L) stop("smoke_n cannot be negative")
smoke_run <- smoke_n > 0L
smoke_selection <- NA_character_
if (smoke_run) {
    ## Sample evenly across the catalog rather than taking its head. The catalog
    ## is sorted by position, so head() always selects the chr1 telomere, where
    ## cis coverage is at its worst: a two-locus smoke there returned two
    ## qc_failed rows, zero complete features, and could never reach Stages
    ## 03-06. Even spacing is deterministic (no RNG, no seed to record),
    ## reproduces exactly for a given smoke_n, spans chromosomes, and still
    ## picks up sparse loci at a realistic rate, so both the complete-feature
    ## and the exclusion paths get exercised.
    smoke_n <- min(smoke_n, nrow(tasks))
    idx <- unique(round(seq(1, nrow(tasks), length.out = smoke_n)))
    tasks <- tasks[idx, , drop = FALSE]
    smoke_selection <- "evenly_spaced_catalog_positions"
}
tasks$task_id <- seq_len(nrow(tasks))
vmrs_per_chunk <- as_int(cli$vmrs_per_chunk, "vmrs_per_chunk")
if (vmrs_per_chunk < 1L) stop("vmrs_per_chunk must be positive")
chunk_manifest <- data.frame(
    chunk_id = ceiling(tasks$task_id / vmrs_per_chunk),
    task_id = tasks$task_id,
    vmr_id = tasks$vmr_id,
    stringsAsFactors = FALSE
)

joint_config <- file.path(module_root, "config", "joint-pve-20260820.tsv")
relative_config <- file.path(repo_root, "config", "local_genetic_control.yml")
threshold_config <- file.path(repo_root, "config", "thresholds.yml")
support_config <- file.path(module_root, "config",
                            "joint-pve-characterized-support.tsv")
for (path in c(joint_config, relative_config, threshold_config,
               support_config)) {
    if (!file.exists(path)) stop("Locked configuration is missing: ", path)
}

dir.create(run_dir, recursive = TRUE)
for (subdir in c("config", "logs", "work", "results/task_rows",
                 "results/combined")) {
    dir.create(file.path(run_dir, subdir), recursive = TRUE)
}
write_tsv(tasks, file.path(run_dir, "config", "task-manifest.tsv"))
write_tsv(chunk_manifest, file.path(run_dir, "config", "chunk-manifest.tsv"))
invisible(file.copy(c(joint_config, relative_config, threshold_config,
                      support_config),
                    file.path(run_dir, "config"), overwrite = FALSE))

git_commit <- tryCatch(
    system2("git", c("-C", repo_root, "rev-parse", "HEAD"), stdout = TRUE)[[1L]],
    error = function(e) NA_character_
)
manifest <- data.frame(
    field = c(
        "run_id", "analysis", "cohort", "region", "started_at",
        "catalog_cohort", "estimation_group", "covar_prefix",
        "upstream_module", "upstream_catalog_run_id",
        "git_commit", "smoke_run", "upstream_vmr_run_id", "vmr_set_id",
        "ordered_donor_checksum", "n_donors", "n_expected_tasks",
        "vmrs_per_chunk", "n_expected_chunks", "smoke_selection",
        "joint_model_run_id", "joint_model_path", "joint_model_sha256",
        "development_features_path", "config_joint_sha256",
        "config_relative_score_sha256", "config_thresholds_sha256",
        "config_characterized_support_sha256"
    ),
    value = c(
        cli$run_id, "02_local_genetic_variance", cohort, region,
        format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        catalog_cohort, estimation_group, covar_prefix,
        upstream_module,
        ## For a cell, the Module 01 run its loci actually came from. Provenance
        ## has to reach past 01b or the locus set is unattributable.
        if (cell$is_estimation_cell) manifest_value("upstream_vmr_catalog")
        else cli$vmr_run_id,
        git_commit,
        toupper(as.character(smoke_run)), cli$vmr_run_id,
        manifest_value("vmr_set_id"), manifest_value("donor_checksum"),
        manifest_value("n_donors"), nrow(tasks), vmrs_per_chunk,
        length(unique(chunk_manifest$chunk_id)), smoke_selection,
        cli$model_run_id,
        normalizePath(model_path), observed_model_sha,
        normalizePath(development_features), sha256_file(joint_config),
        sha256_file(relative_config), sha256_file(threshold_config),
        sha256_file(support_config)
    ),
    stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(run_dir, "manifest.tsv"))
cat(normalizePath(run_dir), "\n", sep = "")
