#!/usr/bin/env Rscript

## Stage 01: calculate the frozen joint estimator's observed-data features for
## one VMR. Every task writes exactly one terminal row, including input-QC and
## computational failures, so Stage 02 can reconcile the full task universe.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
h_dir <- dirname(script_path)
source(file.path(h_dir, "00_functions.R"))
source(file.path(h_dir, "bslmm_pilot_functions.R"))
source(file.path(h_dir, "joint_pve_functions.R"))
source(file.path(h_dir, "observed_locus_io.R"))

cli <- parse_cli(list(
    run_dir = "",
    task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = ""),
    keep_work = "FALSE"
))
if (!nzchar(cli$run_dir) || !nzchar(cli$task_id)) {
    stop("--run-dir and --task-id are required")
}
run_dir <- normalizePath(cli$run_dir)
task_id <- as_int(cli$task_id, "task_id")
keep_work <- as_bool(cli$keep_work, "keep_work")

manifest <- read_tsv(file.path(run_dir, "manifest.tsv"))
mval <- function(field) {
    value <- manifest$value[manifest$field == field]
    if (length(value) != 1L) stop("Run manifest lacks unique field: ", field)
    as.character(value[[1L]])
}
tasks <- read_tsv(file.path(run_dir, "config", "task-manifest.tsv"))
task <- tasks[tasks$task_id == task_id, , drop = FALSE]
if (nrow(task) != 1L) stop("task_id is absent or duplicated: ", task_id)

repo_root <- normalizePath(file.path(h_dir, "..", ".."))
vmr_run_dir <- file.path(
    repo_root, "01_vmr_catalog", "_m", "runs", mval("upstream_vmr_run_id")
)
settings_path <- file.path(run_dir, "config", "joint-pve-20260820.tsv")
settings <- read_joint_settings(settings_path)
threshold_lines <- readLines(file.path(run_dir, "config", "thresholds.yml"),
                             warn = FALSE)
minimum_line <- grep("^[[:space:]]+min_cis_variants:", threshold_lines,
                     value = TRUE)
if (length(minimum_line) != 1L) stop("Cannot resolve min_cis_variants")
min_cis_variants <- as_int(sub(".*:[[:space:]]*", "", minimum_line),
                           "min_cis_variants")

stable_seed <- function(...) {
    text <- paste(..., collapse = "|")
    value <- 104729
    for (byte in utf8ToInt(text)) {
        value <- (value * 131 + byte) %% 2147483629
    }
    as.integer(max(1, value))
}
seed <- stable_seed(mval("run_id"), mval("region"), task$vmr_id,
                    "joint_features")

blank_row <- function() {
    data.frame(
        task_id = task_id,
        cohort = mval("cohort"), region = mval("region"),
        population = mval("cohort"),
        chrom = as.character(task$chrom),
        start = as.integer(task$start), end = as.integer(task$end),
        vmr_id = as.character(task$vmr_id),
        vmr_set_id = as.character(task$vmr_set_id),
        upstream_vmr_run_id = mval("upstream_vmr_run_id"),
        n_cpgs = as.integer(task$n_cpgs),
        n = NA_integer_, samples = NA_integer_,
        num_snps = NA_integer_, n_variants = NA_integer_,
        snps_in_window = NA_integer_, p_eff = NA_real_, ld_metric = NA_real_,
        mean_methylation = NA_real_, methylation_variance = NA_real_,
        bslmm_pve = NA_real_, bslmm_pve_median = NA_real_,
        bslmm_pve_q025 = NA_real_, bslmm_pve_q975 = NA_real_,
        bslmm_h_mean = NA_real_, bslmm_converged = FALSE,
        bslmm_exit_status = NA_integer_, bslmm_elapsed_sec = NA_real_,
        bslmm_n_mcmc = NA_integer_,
        he_h2 = NA_real_, he_se = NA_real_, he_pvalue = NA_real_,
        he_converged = FALSE,
        rho2_oof = NA_real_, r2_oof = NA_real_,
        covariance_ratio_oof = NA_real_,
        score_variance_ratio_oof = NA_real_, en_converged = FALSE,
        feature_complete = FALSE, computational_failure = FALSE,
        terminal_status = NA_character_, exclusion_reason = NA_character_,
        feature_error = NA_character_, feature_seed = seed,
        plink_source = NA_character_, phenotype_source = NA_character_,
        stringsAsFactors = FALSE
    )
}
finish <- function(row, status, reason = NA_character_, error = NA_character_,
                   computational = FALSE) {
    row$terminal_status <- status
    row$exclusion_reason <- reason
    row$feature_error <- error
    row$computational_failure <- computational
    row
}

estimate_task <- function() {
    row <- blank_row()
    locus <- load_observed_locus(
        task = task, cohort = mval("cohort"), vmr_run_dir = vmr_run_dir,
        min_cis_variants = min_cis_variants,
        expected_n = as_int(mval("n_donors"), "n_donors"),
        backing_tag = paste0("lgv-", task_id)
    )
    if (!identical(locus$status, "ok")) {
        if (!is.null(locus$snps_in_window)) {
            row$snps_in_window <- as.integer(locus$snps_in_window)
        }
        return(finish(row, locus$status, locus$reason))
    }
    genotype <- locus$genotype
    y <- locus$y
    covariates <- locus$covariates
    row$snps_in_window <- as.integer(locus$snps_in_window)
    row$num_snps <- row$n_variants <- ncol(genotype)
    row$plink_source <- locus$plink_source
    row$phenotype_source <- locus$phenotype_source
    row$n <- row$samples <- nrow(genotype)
    row$mean_methylation <- mean(y)
    row$methylation_variance <- stats::var(y)

    en <- crossfit_elastic_net(
        genotype = genotype, phenotype = y, covariates = covariates,
        outer_folds = as.integer(settings$outer_folds),
        outer_repeats = as.integer(settings$outer_repeats),
        inner_folds = as.integer(settings$inner_folds),
        alpha_grid = split_numeric(settings$alpha_grid),
        lambda_rule = settings$lambda_rule,
        max_features = as.integer(settings$max_features),
        seed = seed + 17L, keep_predictions = FALSE
    )
    he <- haseman_elston(genotype, y, covariates)
    row$p_eff <- effective_rank_genotype(genotype)
    row$ld_metric <- adjacent_ld_metric(genotype)

    ## GEMMA cannot consume missing BIMBAM dosages. Mean imputation here is a
    ## phenotype-independent full-data operation for the full-data BSLMM
    ## diagnostic; nested EN imputation remains outer-training-only.
    bslmm_genotype <- genotype
    means <- colMeans(bslmm_genotype, na.rm = TRUE)
    for (j in seq_len(ncol(bslmm_genotype))) {
        missing <- is.na(bslmm_genotype[, j])
        if (any(missing)) bslmm_genotype[missing, j] <- means[[j]]
    }
    bslmm_work <- file.path(run_dir, "work", sprintf("vmr-%07d", task_id))
    bslmm <- fit_bslmm_pve(
        genotype = bslmm_genotype,
        phenotype = residualize_phenotype(y, covariates),
        work_dir = bslmm_work,
        gemma_bin = settings$gemma_bin,
        bslmm_mode = as.integer(settings$bslmm_mode),
        burn_in = as.integer(settings$bslmm_burn_in),
        sampling = as.integer(settings$bslmm_sampling),
        rpace = as.integer(settings$bslmm_rpace), seed = seed
    )
    if (!keep_work) unlink(bslmm_work, recursive = TRUE)

    row$bslmm_pve <- bslmm$pve_mean
    row$bslmm_pve_median <- bslmm$pve_median
    row$bslmm_pve_q025 <- bslmm$pve_q025
    row$bslmm_pve_q975 <- bslmm$pve_q975
    row$bslmm_h_mean <- bslmm$h_mean
    row$bslmm_converged <- isTRUE(bslmm$converged)
    row$bslmm_exit_status <- bslmm$exit_status
    row$bslmm_elapsed_sec <- bslmm$elapsed_sec
    row$bslmm_n_mcmc <- bslmm$n_mcmc
    row$he_h2 <- he$he_h2[[1L]]
    row$he_se <- he$he_se[[1L]]
    row$he_pvalue <- he$he_pvalue[[1L]]
    row$he_converged <- isTRUE(he$he_converged[[1L]])
    row$rho2_oof <- en$metrics$rho2_oof[[1L]]
    row$r2_oof <- en$metrics$r2_oof[[1L]]
    row$covariance_ratio_oof <- en$metrics$covariance_ratio_oof[[1L]]
    row$score_variance_ratio_oof <- en$metrics$score_variance_ratio_oof[[1L]]
    row$en_converged <- isTRUE(en$metrics$converged[[1L]])
    required_features <- c(row$bslmm_pve, row$he_h2, row$rho2_oof,
                           row$r2_oof, row$p_eff, row$ld_metric)
    row$feature_complete <- all(is.finite(required_features)) &&
        row$bslmm_converged && row$he_converged && row$en_converged
    if (!row$feature_complete) {
        detail <- paste(na.omit(c(
            if (!row$bslmm_converged) paste0("BSLMM: ", bslmm$error),
            if (!row$he_converged) "HE did not converge",
            if (!row$en_converged) "nested EN did not converge",
            if (any(!is.finite(required_features))) "nonfinite joint feature"
        )), collapse = " | ")
        return(finish(row, "computational_failure", error = detail,
                      computational = TRUE))
    }
    finish(row, "completed")
}

result <- tryCatch(
    estimate_task(),
    error = function(e) finish(
        blank_row(), "computational_failure",
        error = conditionMessage(e), computational = TRUE
    )
)
output <- file.path(run_dir, "results", "task_rows",
                    sprintf("vmr-%07d.tsv", task_id))
write_tsv(result, output)
cat(result$terminal_status, result$vmr_id, "\n")
