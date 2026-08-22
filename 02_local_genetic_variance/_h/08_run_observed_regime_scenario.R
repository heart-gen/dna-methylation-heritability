#!/usr/bin/env Rscript

## Observed-regime grid, stage B: compute the frozen estimator's raw features
## for one chunk of scenarios.
##
## Each scenario keeps a real VMR's genotype, covariates and donor alignment,
## and replaces only the methylation phenotype with one simulated at a known
## true PVE. Every scenario writes exactly one terminal row so stage C can
## reconcile the full scenario universe.

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
    chunk_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = ""),
    keep_work = "FALSE"
))
if (!nzchar(cli$run_dir) || !nzchar(cli$chunk_id)) {
    stop("--run-dir and --chunk-id are required")
}
run_dir <- normalizePath(cli$run_dir)
chunk_id <- as_int(cli$chunk_id, "chunk_id")
keep_work <- as_bool(cli$keep_work, "keep_work")

manifest <- read_tsv(file.path(run_dir, "manifest.tsv"))
mval <- function(field) {
    value <- manifest$value[manifest$field == field]
    if (length(value) != 1L) stop("Run manifest lacks unique field: ", field)
    as.character(value[[1L]])
}
scenarios <- read_tsv(file.path(run_dir, "config", "scenario-manifest.tsv"))
chunks <- read_tsv(file.path(run_dir, "config", "chunk-manifest.tsv"))
wanted <- chunks$scenario_id[chunks$chunk_id == chunk_id]
if (!length(wanted)) stop("No scenarios for chunk ", chunk_id)

repo_root <- normalizePath(file.path(h_dir, "..", ".."))
vmr_run_dir <- file.path(
    repo_root, "01_vmr_catalog", "_m", "runs", mval("upstream_vmr_run_id")
)
threshold_lines <- readLines(file.path(run_dir, "config", "thresholds.yml"),
                             warn = FALSE)
minimum_line <- grep("^[[:space:]]+min_cis_variants:", threshold_lines,
                     value = TRUE)
if (length(minimum_line) != 1L) stop("Cannot resolve min_cis_variants")
min_cis_variants <- as_int(sub(".*:[[:space:]]*", "", minimum_line),
                           "min_cis_variants")
cohort <- mval("cohort")
expected_n <- as_int(mval("n_donors"), "n_donors")

blank_row <- function(scenario) {
    data.frame(
        scenario_id = as.integer(scenario$scenario_id),
        seed = as.numeric(scenario$seed),
        locus_id = as.integer(scenario$locus_id),
        observed_task_id = as.integer(scenario$observed_task_id),
        vmr_id = as.character(scenario$vmr_id),
        chrom = as.character(scenario$chrom),
        start = as.integer(scenario$start), end = as.integer(scenario$end),
        num_snps_stratum = as.integer(scenario$num_snps_stratum),
        p_eff_stratum = as.integer(scenario$p_eff_stratum),
        architecture = as.character(scenario$architecture),
        true_h2 = as.numeric(scenario$true_h2),
        replicate = as.integer(scenario$replicate),
        realized_h2 = NA_real_, num_causal = NA_integer_,
        n = NA_integer_, num_snps = NA_integer_,
        p_eff = NA_real_, ld_metric = NA_real_, mean_maf = NA_real_,
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
        feature_error = NA_character_,
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

## One locus is reused by every scenario in a chunk, so load it once.
locus_cache <- new.env(parent = emptyenv())
get_locus <- function(scenario) {
    key <- as.character(scenario$locus_id)
    if (!is.null(locus_cache[[key]])) return(locus_cache[[key]])
    locus <- load_observed_locus(
        task = scenario, cohort = cohort, vmr_run_dir = vmr_run_dir,
        min_cis_variants = min_cis_variants, expected_n = expected_n,
        backing_tag = paste0("lgvreg-", scenario$locus_id)
    )
    rm(list = ls(locus_cache), envir = locus_cache)
    locus_cache[[key]] <- locus
    locus
}

run_scenario <- function(scenario) {
    row <- blank_row(scenario)
    locus <- get_locus(scenario)
    if (!identical(locus$status, "ok")) {
        return(finish(row, locus$status, locus$reason))
    }
    genotype <- locus$genotype
    covariates <- locus$covariates
    row$n <- nrow(genotype)
    row$num_snps <- ncol(genotype)

    ## The grid is only meaningful if it reproduces the production feature
    ## geometry exactly. Fail closed rather than characterise a different locus.
    if (row$num_snps != as.integer(scenario$observed_num_snps)) {
        return(finish(row, "computational_failure",
                      error = "num_snps differs from the observed run",
                      computational = TRUE))
    }

    set.seed(as.integer(scenario$seed %% 2147483629))
    simulated <- simulate_phenotype_on_observed_genotype(
        genotype = genotype, covariates = covariates,
        h2 = as.numeric(scenario$true_h2),
        architecture = as.character(scenario$architecture)
    )
    y <- simulated$phenotype
    row$realized_h2 <- simulated$realized_h2
    row$num_causal <- length(simulated$causal_index)

    en <- crossfit_elastic_net(
        genotype = genotype, phenotype = y, covariates = covariates,
        outer_folds = as.integer(scenario$outer_folds),
        outer_repeats = as.integer(scenario$outer_repeats),
        inner_folds = as.integer(scenario$inner_folds),
        alpha_grid = split_numeric(scenario$alpha_grid),
        lambda_rule = as.character(scenario$lambda_rule),
        max_features = as.integer(scenario$max_features),
        seed = as.integer(scenario$seed %% 2147483629) + 17L,
        keep_predictions = FALSE
    )
    he <- haseman_elston(genotype, y, covariates)
    row$p_eff <- effective_rank_genotype(genotype)
    row$ld_metric <- adjacent_ld_metric(genotype)
    row$mean_maf <- mean(colMeans(genotype, na.rm = TRUE) / 2)

    bslmm_genotype <- genotype
    means <- colMeans(bslmm_genotype, na.rm = TRUE)
    for (j in seq_len(ncol(bslmm_genotype))) {
        missing <- is.na(bslmm_genotype[, j])
        if (any(missing)) bslmm_genotype[missing, j] <- means[[j]]
    }
    bslmm_work <- file.path(run_dir, "work",
                            sprintf("scenario-%07d", scenario$scenario_id))
    if (dir.exists(bslmm_work)) unlink(bslmm_work, recursive = TRUE)
    bslmm <- fit_bslmm_pve(
        genotype = bslmm_genotype,
        phenotype = residualize_phenotype(y, covariates),
        work_dir = bslmm_work,
        gemma_bin = as.character(scenario$gemma_bin),
        bslmm_mode = as.integer(scenario$bslmm_mode),
        burn_in = as.integer(scenario$bslmm_burn_in),
        sampling = as.integer(scenario$bslmm_sampling),
        rpace = as.integer(scenario$bslmm_rpace),
        seed = as.integer(scenario$seed %% 2147483629)
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

for (scenario_id in wanted) {
    scenario <- scenarios[scenarios$scenario_id == scenario_id, , drop = FALSE]
    if (nrow(scenario) != 1L) stop("scenario_id absent or duplicated: ", scenario_id)
    result <- tryCatch(
        run_scenario(scenario),
        error = function(e) finish(
            blank_row(scenario), "computational_failure",
            error = conditionMessage(e), computational = TRUE
        )
    )
    write_tsv(result, file.path(run_dir, "results", "scenario_rows",
                                sprintf("scenario-%07d.tsv", scenario_id)))
    cat(result$terminal_status, scenario_id, result$vmr_id,
        "h2=", result$true_h2, result$architecture, "\n")
}
