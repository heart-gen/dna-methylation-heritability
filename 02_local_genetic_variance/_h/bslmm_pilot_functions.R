## Helpers for the paired BSLMM vs calibrated elastic-net pilot.
## Estimator-screening only: does not change production Module 02 outputs.

read_pilot_settings <- function(path) {
    table <- read_tsv(path)
    if (!all(c("setting", "value") %in% names(table))) {
        stop("Pilot configuration must contain setting and value columns")
    }
    stats::setNames(as.list(table$value), table$setting)
}

verify_calibration_sha256 <- function(path, expected) {
    expected <- tolower(expected)
    if (!grepl("^[0-9a-f]{64}$", expected)) {
        stop("expected_calibration_sha256 must be a 64-character SHA-256 digest")
    }
    out <- suppressWarnings(tryCatch(
        system2("sha256sum", shQuote(normalizePath(path)),
                stdout = TRUE, stderr = FALSE),
        error = function(e) NA_character_
    ))
    observed <- if (length(out) == 0L || is.na(out[[1L]])) {
        NA_character_
    } else {
        tolower(sub(" .*$", "", out[[1L]]))
    }
    if (is.na(observed)) stop("Could not checksum calibration model: ", path)
    if (!identical(observed, expected)) {
        stop(
            "Calibration model checksum mismatch.\n  expected ", expected,
            "\n  observed ", observed
        )
    }
    invisible(observed)
}

residualize_phenotype <- function(phenotype, covariates = NULL) {
    if (is.null(covariates) || ncol(covariates) == 0L) {
        return(as.numeric(scale(phenotype, scale = FALSE)))
    }
    fit <- stats::lm.fit(
        x = cbind(`(Intercept)` = 1, as.matrix(covariates)),
        y = as.numeric(phenotype)
    )
    as.numeric(fit$residuals)
}

write_bimbam_inputs <- function(outdir, genotype, phenotype, prefix = "locus") {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    geno_path <- file.path(outdir, paste0(prefix, ".geno.txt"))
    pheno_path <- file.path(outdir, paste0(prefix, ".pheno.txt"))
    anno_path <- file.path(outdir, paste0(prefix, ".anno.txt"))
    snp_names <- colnames(genotype)
    if (is.null(snp_names)) snp_names <- paste0("snp_", seq_len(ncol(genotype)))
    geno_lines <- vapply(seq_len(ncol(genotype)), function(j) {
        doses <- sprintf("%.6f", genotype[, j])
        paste(c(snp_names[[j]], "A", "T", doses), collapse = ",")
    }, character(1L))
    writeLines(geno_lines, geno_path)
    writeLines(sprintf("%.10g", as.numeric(phenotype)), pheno_path)
    anno <- data.frame(
        snp = snp_names,
        pos = seq_len(ncol(genotype)),
        chr = 1L,
        stringsAsFactors = FALSE
    )
    utils::write.table(
        anno, anno_path, sep = ",", quote = FALSE, row.names = FALSE,
        col.names = FALSE
    )
    list(geno = geno_path, pheno = pheno_path, anno = anno_path)
}

fit_bslmm_pve <- function(genotype, phenotype, work_dir, gemma_bin,
                          bslmm_mode = 1L, burn_in = 10000L,
                          sampling = 100000L, rpace = 10L, seed = 1L) {
    if (!file.exists(gemma_bin)) stop("GEMMA binary not found: ", gemma_bin)
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    inputs <- write_bimbam_inputs(work_dir, genotype, phenotype, prefix = "locus")
    prefix <- "bslmm"
    args <- c(
        "-g", inputs$geno,
        "-p", inputs$pheno,
        "-a", inputs$anno,
        "-bslmm", as.character(bslmm_mode),
        "-w", as.character(burn_in),
        "-s", as.character(sampling),
        "-rpace", as.character(rpace),
        "-seed", as.character(seed),
        "-outdir", work_dir,
        "-o", prefix
    )
    log_path <- file.path(work_dir, "gemma.log")
    started <- proc.time()[["elapsed"]]
    status <- system2(
        gemma_bin, args = args, stdout = log_path, stderr = log_path
    )
    elapsed <- proc.time()[["elapsed"]] - started
    hyp_path <- file.path(work_dir, paste0(prefix, ".hyp.txt"))
    if (!identical(as.integer(status), 0L) || !file.exists(hyp_path)) {
        return(list(
            converged = FALSE,
            pve_mean = NA_real_,
            pve_median = NA_real_,
            pve_q025 = NA_real_,
            pve_q975 = NA_real_,
            h_mean = NA_real_,
            n_mcmc = 0L,
            exit_status = as.integer(status),
            elapsed_sec = elapsed,
            hyp_path = hyp_path,
            error = if (file.exists(log_path)) {
                paste(readLines(log_path, warn = FALSE), collapse = " | ")
            } else {
                "GEMMA failed without a log"
            }
        ))
    }
    hyp <- utils::read.table(
        hyp_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
    )
    names(hyp) <- tolower(gsub("[^A-Za-z0-9_]+", "", names(hyp)))
    if (!"pve" %in% names(hyp)) {
        return(list(
            converged = FALSE,
            pve_mean = NA_real_,
            pve_median = NA_real_,
            pve_q025 = NA_real_,
            pve_q975 = NA_real_,
            h_mean = NA_real_,
            n_mcmc = nrow(hyp),
            exit_status = as.integer(status),
            elapsed_sec = elapsed,
            hyp_path = hyp_path,
            error = "hyp.txt lacks a pve column"
        ))
    }
    pve <- as.numeric(hyp$pve)
    pve <- pve[is.finite(pve)]
    h_col <- if ("h" %in% names(hyp)) as.numeric(hyp$h) else rep(NA_real_, nrow(hyp))
    h_col <- h_col[is.finite(h_col)]
    list(
        converged = length(pve) > 0L,
        pve_mean = if (length(pve)) mean(pve) else NA_real_,
        pve_median = if (length(pve)) stats::median(pve) else NA_real_,
        pve_q025 = if (length(pve)) {
            as.numeric(stats::quantile(pve, 0.025, names = FALSE))
        } else {
            NA_real_
        },
        pve_q975 = if (length(pve)) {
            as.numeric(stats::quantile(pve, 0.975, names = FALSE))
        } else {
            NA_real_
        },
        h_mean = if (length(h_col)) mean(h_col) else NA_real_,
        n_mcmc = length(pve),
        exit_status = as.integer(status),
        elapsed_sec = elapsed,
        hyp_path = hyp_path,
        error = NA_character_
    )
}

run_paired_en_bslmm <- function(scenario, calibration_model, gemma_bin,
                                bslmm_mode, burn_in, sampling, rpace,
                                work_dir, keep_work = FALSE) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(scenario$seed)
    simulated <- simulate_locus(
        n = scenario$n,
        p = scenario$num_snps,
        ld_rho = scenario$ld_rho,
        h2 = scenario$true_h2,
        architecture = scenario$architecture
    )

    en_error <- NA_character_
    fit <- tryCatch(
        crossfit_elastic_net(
            genotype = simulated$genotype,
            phenotype = simulated$phenotype,
            covariates = simulated$covariates,
            outer_folds = scenario$outer_folds,
            outer_repeats = scenario$outer_repeats,
            inner_folds = scenario$inner_folds,
            alpha_grid = split_numeric(scenario$alpha_grid),
            lambda_rule = scenario$lambda_rule,
            max_features = scenario$max_features,
            seed = scenario$seed + 17L,
            keep_predictions = FALSE
        ),
        error = function(e) {
            en_error <<- conditionMessage(e)
            NULL
        }
    )
    if (is.null(fit)) {
        en_metrics <- data.frame(
            n = scenario$n,
            num_snps = scenario$num_snps,
            outer_folds = scenario$outer_folds,
            outer_repeats = scenario$outer_repeats,
            inner_folds = scenario$inner_folds,
            max_features = scenario$max_features,
            r2_oof = NA_real_,
            rho2_oof = NA_real_,
            covariance_ratio_oof = NA_real_,
            score_variance_ratio_oof = NA_real_,
            calibration_slope_oof = NA_real_,
            mean_fold_score_variance_ratio = NA_real_,
            mean_nonzero_snps = NA_real_,
            converged = FALSE
        )
    } else {
        en_metrics <- fit$metrics
    }
    en_metrics$ld_metric <- simulated$ld_metric
    he_metrics <- haseman_elston(
        simulated$genotype, simulated$phenotype, simulated$covariates
    )
    calibrated <- tryCatch(
        predict_calibration(calibration_model, cbind(en_metrics, he_metrics)),
        error = function(e) {
            en_error <<- paste(
                na.omit(c(en_error, conditionMessage(e))), collapse = " | "
            )
            data.frame(
                calibration_stratum = NA_character_,
                calibration_distance = NA_real_,
                h2_en_forward = NA_real_,
                h2_he_weight = NA_real_,
                h2_en_calibrated_unbounded = NA_real_,
                h2_en_calibrated = NA_real_,
                h2_calibration_lower = NA_real_,
                h2_calibration_upper = NA_real_,
                null_raw_threshold_95 = NA_real_,
                calibration_upper_bound = NA_real_,
                h2_upper_boundary_hit = NA,
                positive_signal = NA,
                calibration_status = "raw_metric_unavailable",
                stringsAsFactors = FALSE
            )
        }
    )

    residual_y <- residualize_phenotype(
        simulated$phenotype, simulated$covariates
    )
    bslmm_dir <- file.path(work_dir, sprintf("scenario-%07d", scenario$scenario_id))
    if (dir.exists(bslmm_dir)) unlink(bslmm_dir, recursive = TRUE)
    bslmm <- fit_bslmm_pve(
        genotype = simulated$genotype,
        phenotype = residual_y,
        work_dir = bslmm_dir,
        gemma_bin = gemma_bin,
        bslmm_mode = bslmm_mode,
        burn_in = burn_in,
        sampling = sampling,
        rpace = rpace,
        seed = scenario$seed
    )
    if (!isTRUE(keep_work)) {
        unlink(bslmm_dir, recursive = TRUE)
    }

    truth <- scenario$true_h2
    en_est <- calibrated$h2_en_calibrated[[1L]]
    bslmm_est <- bslmm$pve_mean
    data.frame(
        scenario_id = scenario$scenario_id,
        seed = scenario$seed,
        n = scenario$n,
        num_snps = scenario$num_snps,
        ld_rho = scenario$ld_rho,
        architecture = scenario$architecture,
        true_h2 = truth,
        replicate = scenario$replicate,
        realized_h2 = simulated$realized_h2,
        num_causal = length(simulated$causal_index),
        ld_metric = simulated$ld_metric,
        mean_maf = mean(colMeans(simulated$genotype) / 2),
        en_converged = isTRUE(en_metrics$converged[[1L]]),
        en_error = en_error,
        r2_oof = en_metrics$r2_oof[[1L]],
        rho2_oof = en_metrics$rho2_oof[[1L]],
        he_h2 = he_metrics$he_h2[[1L]],
        h2_en_calibrated = en_est,
        h2_en_calibrated_unbounded = calibrated$h2_en_calibrated_unbounded[[1L]],
        h2_calibration_lower = calibrated$h2_calibration_lower[[1L]],
        h2_calibration_upper = calibrated$h2_calibration_upper[[1L]],
        h2_upper_boundary_hit = calibrated$h2_upper_boundary_hit[[1L]],
        calibration_status = calibrated$calibration_status[[1L]],
        en_error_vs_truth = en_est - truth,
        en_covered = is.finite(en_est) &&
            is.finite(calibrated$h2_calibration_lower[[1L]]) &&
            is.finite(calibrated$h2_calibration_upper[[1L]]) &&
            truth >= calibrated$h2_calibration_lower[[1L]] &&
            truth <= calibrated$h2_calibration_upper[[1L]],
        en_failed = !(isTRUE(en_metrics$converged[[1L]]) && is.finite(en_est)),
        en_boundary_or_extrap = isTRUE(calibrated$h2_upper_boundary_hit[[1L]]) ||
            identical(calibrated$calibration_status[[1L]], "outside_design_domain") ||
            identical(calibrated$calibration_status[[1L]], "raw_metric_extrapolation"),
        bslmm_converged = isTRUE(bslmm$converged),
        bslmm_exit_status = bslmm$exit_status,
        bslmm_error = bslmm$error,
        bslmm_elapsed_sec = bslmm$elapsed_sec,
        bslmm_n_mcmc = bslmm$n_mcmc,
        bslmm_pve = bslmm_est,
        bslmm_pve_median = bslmm$pve_median,
        bslmm_pve_q025 = bslmm$pve_q025,
        bslmm_pve_q975 = bslmm$pve_q975,
        bslmm_h_mean = bslmm$h_mean,
        bslmm_error_vs_truth = bslmm_est - truth,
        bslmm_covered = is.finite(bslmm_est) &&
            is.finite(bslmm$pve_q025) && is.finite(bslmm$pve_q975) &&
            truth >= bslmm$pve_q025 && truth <= bslmm$pve_q975,
        bslmm_failed = !isTRUE(bslmm$converged) || !is.finite(bslmm_est),
        stringsAsFactors = FALSE
    )
}

run_bslmm_only <- function(scenario, gemma_bin, bslmm_mode, burn_in,
                           sampling, rpace, work_dir, keep_work = FALSE,
                           positive_signal_rule = "mcmc_lower_gt_0") {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(scenario$seed)
    simulated <- simulate_locus(
        n = scenario$n,
        p = scenario$num_snps,
        ld_rho = scenario$ld_rho,
        h2 = scenario$true_h2,
        architecture = scenario$architecture
    )
    residual_y <- residualize_phenotype(
        simulated$phenotype, simulated$covariates
    )
    bslmm_dir <- file.path(work_dir, sprintf("scenario-%07d", scenario$scenario_id))
    if (dir.exists(bslmm_dir)) unlink(bslmm_dir, recursive = TRUE)
    bslmm <- fit_bslmm_pve(
        genotype = simulated$genotype,
        phenotype = residual_y,
        work_dir = bslmm_dir,
        gemma_bin = gemma_bin,
        bslmm_mode = bslmm_mode,
        burn_in = burn_in,
        sampling = sampling,
        rpace = rpace,
        seed = scenario$seed
    )
    if (!isTRUE(keep_work)) unlink(bslmm_dir, recursive = TRUE)

    truth <- scenario$true_h2
    est <- bslmm$pve_mean
    failed <- !isTRUE(bslmm$converged) || !is.finite(est)
    in_unit <- is.finite(est) && est >= 0 && est <= 1
    ## BSLMM has no EN-style design-distance domain. Operational domain for the
    ## Module 02 analog gate: converged, finite, and in [0, 1].
    within_domain <- !failed && in_unit
    positive_signal <- if (identical(positive_signal_rule, "mcmc_lower_gt_0")) {
        is.finite(bslmm$pve_q025) && bslmm$pve_q025 > 0
    } else {
        stop("Unknown positive_signal_rule: ", positive_signal_rule)
    }
    covered <- is.finite(est) && is.finite(bslmm$pve_q025) &&
        is.finite(bslmm$pve_q975) &&
        truth >= bslmm$pve_q025 && truth <= bslmm$pve_q975

    data.frame(
        scenario_id = scenario$scenario_id,
        split = if ("split" %in% names(scenario)) scenario$split[[1L]] else "evaluation",
        seed = scenario$seed,
        n = scenario$n,
        num_snps = scenario$num_snps,
        ld_rho = scenario$ld_rho,
        architecture = scenario$architecture,
        true_h2 = truth,
        replicate = scenario$replicate,
        realized_h2 = simulated$realized_h2,
        num_causal = length(simulated$causal_index),
        ld_metric = simulated$ld_metric,
        mean_maf = mean(colMeans(simulated$genotype) / 2),
        bslmm_converged = isTRUE(bslmm$converged),
        bslmm_exit_status = bslmm$exit_status,
        bslmm_error = bslmm$error,
        bslmm_elapsed_sec = bslmm$elapsed_sec,
        bslmm_n_mcmc = bslmm$n_mcmc,
        bslmm_pve = est,
        bslmm_pve_median = bslmm$pve_median,
        bslmm_pve_q025 = bslmm$pve_q025,
        bslmm_pve_q975 = bslmm$pve_q975,
        bslmm_h_mean = bslmm$h_mean,
        error = est - truth,
        covered = covered,
        positive_signal = positive_signal,
        failed = failed,
        within_domain = within_domain,
        calibration_status = if (within_domain) {
            "within_domain"
        } else if (failed) {
            "computational_failure"
        } else {
            "outside_unit_interval"
        },
        stringsAsFactors = FALSE
    )
}
