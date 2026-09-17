#!/usr/bin/env Rscript

## Stage 15: open a RESCORE run from a sealed, accepted observed production run.
##
## Why this exists. The characterized-support table gates eligibility only: it
## sets joint_pve_domain_status in Stage 03, which Stage 04 turns into
## `local_genetic_control_eligible` and hence into the midrank denominator. It
## never enters the estimator. pve_cis_joint_unbounded is produced by the
## checksum-pinned frozen calibrator from observed-joint-features.tsv alone, so
## a support-table correction changes WHICH loci are scored and the percentile
## they receive, and changes no estimate.
##
## That distinction is what this stage exploits. The six 20260823 arm runs were
## sealed against a support table pooled across cells (sha f52026944b...), so
## the AA arm inherited a p_eff floor of 2.058 from `all_individuals` instead of
## its own 7.079, and scored 82/83/91 loci (0.72-0.98%) that AA's own regime
## grids never characterised. Re-estimating those runs would burn ~1,900 array
## tasks each to reproduce byte-identical features. Instead this stage opens a
## new immutable run that CARRIES the sealed features forward and re-runs only
## Stages 03-06 against the current per-cell support table.
##
## Runs stay immutable (AGENTS.md 5.2): the source run is read-only and is not
## touched. The rescore is a new run with its own ID, its own manifest, and its
## own Stage 05 gate, and like any other run it is not accepted until a human
## records it in the module README (AGENTS.md 6).

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
h_dir <- dirname(script_path)
source(file.path(h_dir, "00_functions.R"))

repo_root <- normalizePath(file.path(h_dir, "..", ".."))
module_root <- file.path(repo_root, "02_local_genetic_variance")
cli <- parse_cli(list(run_id = "", source_run_id = "", rescore_reason = "",
                      runs_root = ""))
required_cli <- c("run_id", "source_run_id", "rescore_reason")
missing_cli <- required_cli[!vapply(cli[required_cli], nzchar, logical(1L))]
if (length(missing_cli)) {
    stop("Required argument(s): ", paste(missing_cli, collapse = ", "))
}
runs_root <- if (nzchar(cli$runs_root)) cli$runs_root else
    file.path(module_root, "_m", "runs")
source_dir <- file.path(runs_root, cli$source_run_id)
run_dir <- file.path(runs_root, cli$run_id)
if (!dir.exists(source_dir)) stop("Source run does not exist: ", source_dir)
if (dir.exists(run_dir)) stop("Rescore run directory already exists: ", run_dir)
if (identical(cli$run_id, cli$source_run_id)) {
    stop("A rescore run must have its own run ID")
}

src <- read_tsv(file.path(source_dir, "manifest.tsv"))
sval <- function(field, required = TRUE) {
    value <- src$value[src$field == field]
    if (length(value) != 1L) {
        if (!required) return(NA_character_)
        stop("Source manifest lacks unique field: ", field)
    }
    as.character(value[[1L]])
}

## The source must be a finished, passing, non-smoke production run. Rescoring
## an open run would race Stage 02; rescoring a failed or smoke run would
## launder it into something that looks acceptable.
if (is.na(sval("finished_at", FALSE))) {
    stop("Source run is not sealed (no finished_at): ", cli$source_run_id)
}
if (!identical(sval("observed_score_decision"),
               "PASS_RELATIVE_SCORE_OBSERVED_QC")) {
    stop("Source run did not pass its own Stage 05 gate: ",
         sval("observed_score_decision"))
}
if (!identical(toupper(sval("smoke_run")), "FALSE")) {
    stop("Refusing to rescore a smoke run: ", cli$source_run_id)
}
## And it must be recorded as accepted by a human, on the same table Module 02
## already publishes. A rescore inherits the source's scientific standing, so
## the source has to have some.
## Module 02's acceptance table writes the run ID bare; Module 01 and 01b wrap
## it in backticks. Match the ID as a table cell under either convention rather
## than assuming one, so this check cannot pass or fail on formatting.
readme <- readLines(file.path(module_root, "README.md"), warn = FALSE)
accepted <- readme[grepl(paste0("|", cli$source_run_id, "|"),
                         gsub("[` ]", "", readme), fixed = TRUE)]
if (length(accepted) != 1L ||
    !grepl("PASS_RELATIVE_SCORE_OBSERVED_QC", accepted, fixed = TRUE)) {
    stop("Source run is not recorded as accepted in the Module 02 README: ",
         cli$source_run_id)
}

features_src <- file.path(source_dir, "results", "combined",
                          "observed-joint-features.tsv")
## Stage 05 reconciles task completion from a table Stage 02 writes. A rescore
## skips Stage 02, but the reconciliation is a fact about the source run's
## array, and the features carried forward are exactly the ones it describes,
## so it is carried alongside them rather than recomputed or waived.
recon_src <- file.path(source_dir, "results", "combined",
                       "task-reconciliation.tsv")
for (path in c(features_src, recon_src)) {
    if (!file.exists(path)) stop("Source run lacks a carried input: ", path)
}
sha256_file <- function(path) {
    out <- system2("sha256sum", normalizePath(path), stdout = TRUE)
    if (!length(out)) stop("Could not checksum ", path)
    tolower(sub(" .*$", "", out[[1L]]))
}
features_sha <- sha256_file(features_src)

## Current locked configuration -- the whole point is that the support table is
## now the corrected per-cell one, so it is read fresh rather than inherited.
joint_config <- file.path(module_root, "config", "joint-pve-20260820.tsv")
relative_config <- file.path(repo_root, "config", "local_genetic_control.yml")
threshold_config <- file.path(repo_root, "config", "thresholds.yml")
support_config <- file.path(module_root, "config",
                            "joint-pve-characterized-support.tsv")
for (path in c(joint_config, relative_config, threshold_config,
               support_config)) {
    if (!file.exists(path)) stop("Locked configuration is missing: ", path)
}
support <- read_tsv(support_config)
if (!"cell" %in% names(support)) {
    stop("Support table is not per-cell; rescoring against a pooled table ",
         "would reinstate the cross-cell borrowing this stage exists to undo")
}
run_cell <- sval("cohort")
if (!any(as.character(support$cell) == run_cell)) {
    stop("Per-cell support table has no rows for cell: ", run_cell)
}
support_sha <- sha256_file(support_config)
if (identical(support_sha, sval("config_characterized_support_sha256"))) {
    stop("Support table is unchanged from the source run; a rescore would ",
         "reproduce it exactly. Nothing to do.")
}

dir.create(run_dir, recursive = TRUE)
for (subdir in c("config", "logs", "work", "results/task_rows",
                 "results/combined")) {
    dir.create(file.path(run_dir, subdir), recursive = TRUE)
}
## Task and chunk manifests come from the source: the locus set is fixed by
## definition, and Stage 05's reconciliation counts against it.
for (f in c("task-manifest.tsv", "chunk-manifest.tsv")) {
    ok <- file.copy(file.path(source_dir, "config", f),
                    file.path(run_dir, "config", f), overwrite = FALSE)
    if (!ok) stop("Could not carry forward ", f)
}
invisible(file.copy(c(joint_config, relative_config, threshold_config,
                      support_config),
                    file.path(run_dir, "config"), overwrite = FALSE))
## The features are carried, not recomputed. This is the reuse that makes a
## rescore cheap, and the checksum is recorded so the reuse is auditable.
for (path in c(features_src, recon_src)) {
    if (!file.copy(path, file.path(run_dir, "results", "combined",
                                   basename(path)), overwrite = FALSE)) {
        stop("Could not carry forward ", basename(path))
    }
}

git_commit <- tryCatch(
    system2("git", c("-C", repo_root, "rev-parse", "HEAD"), stdout = TRUE)[[1L]],
    error = function(e) NA_character_
)
carry <- c("analysis", "cohort", "region", "catalog_cohort",
           "estimation_group", "covar_prefix", "upstream_module",
           "upstream_catalog_run_id", "smoke_run", "upstream_vmr_run_id",
           "vmr_set_id", "ordered_donor_checksum", "n_donors",
           "n_expected_tasks", "vmrs_per_chunk", "n_expected_chunks",
           "smoke_selection", "joint_model_run_id", "joint_model_path",
           "joint_model_sha256", "development_features_path",
           "config_joint_sha256", "config_relative_score_sha256",
           "config_thresholds_sha256")
carried <- vapply(carry, function(f) sval(f, required = FALSE), character(1L))
manifest <- data.frame(
    field = c("run_id", "run_kind", "rescore_of",
              "rescore_reason", "rescore_source_decision",
              "rescore_source_support_sha256",
              "carried_joint_features_sha256", "started_at", "git_commit",
              carry, "config_characterized_support_sha256"),
    value = c(cli$run_id, "observed_rescore", cli$source_run_id,
              cli$rescore_reason, sval("observed_score_decision"),
              sval("config_characterized_support_sha256"),
              features_sha,
              format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), git_commit,
              unname(carried), support_sha),
    stringsAsFactors = FALSE
)
if (anyDuplicated(manifest$field)) {
    stop("Duplicate manifest field in rescore run")
}
write_tsv(manifest, file.path(run_dir, "manifest.tsv"))
cat("Opened rescore run", cli$run_id, "from", cli$source_run_id, "\n")
cat("  cell:", run_cell, " region:", sval("region"), "\n")
cat("  carried features sha256:", features_sha, "\n")
cat("  support sha256:", sval("config_characterized_support_sha256"),
    "->", support_sha, "\n")
cat("  next: Stages 03-06 (no array; features are carried)\n")
