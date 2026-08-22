#!/usr/bin/env Rscript

## Observed-regime grid, stage A: open the run and freeze its scenario universe.
##
## The frozen joint model was calibrated on AR(1) simulated genotypes whose
## effective dimension (p_eff) is far higher than any real cis-window at
## n = 153. This grid removes that mismatch by simulating only the phenotype
## and keeping the real genotype, covariates and donor alignment, so n,
## num_snps, LD and p_eff are the observed values by construction.
##
## Nothing here refits, retunes or reselects any estimator.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
h_dir <- dirname(script_path)
source(file.path(h_dir, "00_functions.R"))
source(file.path(h_dir, "joint_pve_functions.R"))

repo_root <- normalizePath(file.path(h_dir, "..", ".."))
module_root <- file.path(repo_root, "02_local_genetic_variance")
cli <- parse_cli(list(
    config = file.path(module_root, "config", "observed-regime-20260822.tsv"),
    run_id = "",
    runs_root = ""
))
settings <- read_joint_settings(cli$config)
sval <- function(field) {
    value <- settings[[field]]
    if (is.null(value) || !nzchar(as.character(value))) {
        stop("Observed-regime config lacks field: ", field)
    }
    as.character(value)
}
run_id <- if (nzchar(cli$run_id)) cli$run_id else sval("run_id")
if (!grepl("^lgv-observed-regime-[0-9]{8}[a-z]?$", run_id)) {
    stop("run_id must match ^lgv-observed-regime-[0-9]{8}[a-z]?$")
}
if (!identical(sval("model_changed"), "FALSE") ||
    !identical(sval("criteria_changed"), "FALSE")) {
    stop("Observed-regime grid may not change the frozen model or its gates")
}

runs_root <- if (nzchar(cli$runs_root)) cli$runs_root else
    file.path(module_root, "_m", "runs")
run_dir <- file.path(runs_root, run_id)
if (file.exists(run_dir)) stop("Run directory already exists: ", run_dir)

observed_run_dir <- file.path(runs_root, sval("observed_run_id"))
observed_manifest_path <- file.path(observed_run_dir, "manifest.tsv")
observed_features_path <- file.path(
    observed_run_dir, "results", "combined", "observed-joint-features.tsv"
)
for (path in c(observed_manifest_path, observed_features_path)) {
    if (!file.exists(path)) stop("Missing observed-run input: ", path)
}
observed_manifest <- read_tsv(observed_manifest_path)
oval <- function(field) {
    value <- observed_manifest$value[observed_manifest$field == field]
    if (length(value) != 1L) stop("Observed manifest lacks unique field: ", field)
    as.character(value[[1L]])
}

features <- read_tsv(observed_features_path)
eligible <- features[features$terminal_status %in% "completed" &
                     features$feature_complete %in% TRUE, , drop = FALSE]
if (!nrow(eligible)) stop("Observed run contributed no complete feature rows")
for (field in c("task_id", "vmr_id", "chrom", "start", "end", "n_cpgs",
                "num_snps", "p_eff", "ld_metric", "n")) {
    if (!field %in% names(eligible)) {
        stop("Observed features lack field: ", field)
    }
}

## Stratify the real locus universe on the two features that drive the
## mismatch, so the grid spans the whole observed regime rather than its mode.
strata_bins <- function(x, k) {
    breaks <- stats::quantile(as.numeric(x), probs = seq(0, 1, length.out = k + 1L),
                              na.rm = TRUE, type = 7)
    breaks[[1L]] <- -Inf
    breaks[[length(breaks)]] <- Inf
    if (anyDuplicated(breaks)) stop("Stratification breaks are not unique")
    as.integer(cut(as.numeric(x), breaks = breaks, labels = FALSE,
                   include.lowest = TRUE))
}
k_snps <- as_int(sval("num_snps_strata"), "num_snps_strata")
k_peff <- as_int(sval("p_eff_strata"), "p_eff_strata")
per_stratum <- as_int(sval("loci_per_stratum"), "loci_per_stratum")
eligible$num_snps_stratum <- strata_bins(eligible$num_snps, k_snps)
eligible$p_eff_stratum <- strata_bins(eligible$p_eff, k_peff)

base_seed <- as_int(sval("base_seed"), "base_seed") +
    as_int(sval("regime_seed_offset"), "regime_seed_offset")
set.seed(base_seed)
locus_records <- list()
for (a in seq_len(k_snps)) {
    for (b in seq_len(k_peff)) {
        pool <- eligible[eligible$num_snps_stratum == a &
                         eligible$p_eff_stratum == b, , drop = FALSE]
        if (!nrow(pool)) next
        take <- min(per_stratum, nrow(pool))
        pick <- sort(sample(seq_len(nrow(pool)), take, replace = FALSE))
        locus_records[[length(locus_records) + 1L]] <- pool[pick, , drop = FALSE]
    }
}
loci <- do.call(rbind, locus_records)
if (!nrow(loci)) stop("Locus stratification selected nothing")
loci <- loci[order(loci$task_id), , drop = FALSE]
loci$locus_id <- seq_len(nrow(loci))
locus_manifest <- data.frame(
    locus_id = loci$locus_id,
    observed_task_id = as.integer(loci$task_id),
    vmr_id = as.character(loci$vmr_id),
    chrom = as.character(loci$chrom),
    start = as.integer(loci$start),
    end = as.integer(loci$end),
    n_cpgs = as.integer(loci$n_cpgs),
    num_snps_stratum = as.integer(loci$num_snps_stratum),
    p_eff_stratum = as.integer(loci$p_eff_stratum),
    observed_n = as.integer(loci$n),
    observed_num_snps = as.integer(loci$num_snps),
    observed_p_eff = as.numeric(loci$p_eff),
    observed_ld_metric = as.numeric(loci$ld_metric),
    stringsAsFactors = FALSE
)

architectures <- split_character(sval("architectures"))
h2_values <- split_numeric(sval("h2_values"))
replicates <- as_int(sval("replicates_per_cell"), "replicates_per_cell")
grid <- expand.grid(
    replicate = seq_len(replicates),
    true_h2 = h2_values,
    architecture = architectures,
    locus_id = locus_manifest$locus_id,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
grid <- grid[order(grid$locus_id, grid$architecture, grid$true_h2,
                   grid$replicate), , drop = FALSE]
scenarios <- merge(grid, locus_manifest, by = "locus_id", sort = FALSE)
scenarios <- scenarios[order(scenarios$locus_id, scenarios$architecture,
                             scenarios$true_h2, scenarios$replicate), ,
                       drop = FALSE]
scenarios$scenario_id <- seq_len(nrow(scenarios))
scenarios$seed <- base_seed + scenarios$scenario_id * 1009L
if (anyDuplicated(scenarios$scenario_id)) stop("Duplicate scenario_id")
if (anyDuplicated(scenarios$seed)) stop("Duplicate scenario seed")

for (field in c("outer_folds", "outer_repeats", "inner_folds", "max_features",
                "bslmm_mode", "bslmm_burn_in", "bslmm_sampling",
                "bslmm_rpace")) {
    scenarios[[field]] <- as_int(sval(field), field)
}
for (field in c("alpha_grid", "lambda_rule", "gemma_bin")) {
    scenarios[[field]] <- sval(field)
}
scenarios <- scenarios[, c(
    "scenario_id", "seed", "locus_id", "observed_task_id", "vmr_id", "chrom",
    "start", "end", "n_cpgs", "num_snps_stratum", "p_eff_stratum",
    "observed_n", "observed_num_snps", "observed_p_eff", "observed_ld_metric",
    "architecture", "true_h2", "replicate", "outer_folds", "outer_repeats",
    "inner_folds", "alpha_grid", "lambda_rule", "max_features", "gemma_bin",
    "bslmm_mode", "bslmm_burn_in", "bslmm_sampling", "bslmm_rpace"
), drop = FALSE]

per_chunk <- as_int(sval("scenarios_per_chunk"), "scenarios_per_chunk")
chunk_manifest <- data.frame(
    chunk_id = ((seq_len(nrow(scenarios)) - 1L) %/% per_chunk) + 1L,
    scenario_id = scenarios$scenario_id,
    stringsAsFactors = FALSE
)

dir.create(run_dir, recursive = TRUE)
for (subdir in c("config", "logs", "work", "results/scenario_rows",
                 "results/combined")) {
    dir.create(file.path(run_dir, subdir), recursive = TRUE)
}
threshold_config <- file.path(repo_root, "config", "thresholds.yml")
invisible(file.copy(c(normalizePath(cli$config), threshold_config),
                    file.path(run_dir, "config"), overwrite = FALSE))
write_tsv(locus_manifest, file.path(run_dir, "config", "locus-manifest.tsv"))
write_tsv(scenarios, file.path(run_dir, "config", "scenario-manifest.tsv"))
write_tsv(chunk_manifest, file.path(run_dir, "config", "chunk-manifest.tsv"))

sha256 <- function(path) {
    tolower(sub(" .*$", "", system2("sha256sum", normalizePath(path),
                                    stdout = TRUE)[[1L]]))
}
model_path <- oval("joint_model_path")
model_sha <- sha256(model_path)
if (!identical(model_sha, tolower(oval("joint_model_sha256")))) {
    stop("Frozen joint-model checksum differs from the observed run")
}
manifest <- data.frame(
    field = c(
        "run_id", "analysis", "grid_kind", "cohort", "region", "started_at",
        "git_commit", "observed_run_id", "upstream_vmr_run_id", "vmr_set_id",
        "ordered_donor_checksum", "n_donors", "n_loci", "n_scenarios",
        "scenarios_per_chunk", "n_expected_chunks", "joint_model_run_id",
        "joint_model_path", "joint_model_sha256",
        "development_features_path", "config_observed_regime_sha256",
        "base_seed", "regime_seed_offset", "model_changed", "criteria_changed"
    ),
    value = c(
        run_id, "02_local_genetic_variance", "observed_regime_real_genotype",
        oval("cohort"), oval("region"),
        format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        system2("git", c("-C", repo_root, "rev-parse", "HEAD"), stdout = TRUE),
        sval("observed_run_id"), oval("upstream_vmr_run_id"),
        oval("vmr_set_id"), oval("ordered_donor_checksum"), oval("n_donors"),
        nrow(locus_manifest), nrow(scenarios), per_chunk,
        max(chunk_manifest$chunk_id), oval("joint_model_run_id"),
        model_path, model_sha, oval("development_features_path"),
        sha256(cli$config), sval("base_seed"), sval("regime_seed_offset"),
        "FALSE", "FALSE"
    ),
    stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(run_dir, "manifest.tsv"))
cat(normalizePath(run_dir), "\n")
