## Helpers for the final prespecified joint BSLMM/EN/HE PVE experiment.
## The estimator family is locked in config/FINAL_JOINT_PVE_STRATEGY.md.

read_joint_settings <- function(path) {
    tab <- read_tsv(path)
    if (!all(c("setting", "value") %in% names(tab))) {
        stop("Joint PVE config must contain setting and value columns")
    }
    stats::setNames(as.list(tab$value), tab$setting)
}

split_character <- function(value) {
    trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1L]])
}

effective_rank_genotype <- function(genotype) {
    genotype <- as.matrix(genotype)
    storage.mode(genotype) <- "double"
    if (anyNA(genotype)) {
        means <- colMeans(genotype, na.rm = TRUE)
        means[!is.finite(means)] <- 0
        for (j in seq_len(ncol(genotype))) {
            miss <- is.na(genotype[, j])
            if (any(miss)) genotype[miss, j] <- means[[j]]
        }
    }
    variances <- apply(genotype, 2L, stats::var)
    keep <- is.finite(variances) & variances > 1e-8
    if (!any(keep)) return(NA_real_)
    z <- scale(genotype[, keep, drop = FALSE])
    k <- tcrossprod(z) / ncol(z)
    trace_k <- sum(diag(k))
    trace_k2 <- sum(k * k)
    if (!is.finite(trace_k) || !is.finite(trace_k2) || trace_k2 <= 0) {
        return(NA_real_)
    }
    as.numeric(trace_k^2 / trace_k2)
}

joint_missing_en_metrics <- function(scenario) {
    data.frame(
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
}

assert_source_identity <- function(source, scenario, simulated) {
    fields <- c("scenario_id", "seed", "n", "num_snps", "architecture",
                "true_h2", "replicate")
    for (field in fields) {
        observed <- as.character(source[[field]][[1L]])
        expected <- as.character(scenario[[field]][[1L]])
        if (!identical(observed, expected)) {
            stop("Development BSLMM source mismatch for ", field,
                 ": expected ", expected, " observed ", observed)
        }
    }
    numeric_checks <- c(ld_rho = scenario$ld_rho[[1L]],
                        ld_metric = simulated$ld_metric)
    for (field in names(numeric_checks)) {
        if (!isTRUE(all.equal(as.numeric(source[[field]][[1L]]),
                              as.numeric(numeric_checks[[field]]),
                              tolerance = 1e-10))) {
            stop("Development BSLMM source mismatch for ", field)
        }
    }
    invisible(TRUE)
}

run_joint_pve_features <- function(scenario, development_bslmm_root = "",
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

    errors <- character()
    en_fit <- tryCatch(
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
            errors <<- c(errors, paste0("EN: ", conditionMessage(e)))
            NULL
        }
    )
    en <- if (is.null(en_fit)) joint_missing_en_metrics(scenario) else en_fit$metrics
    he <- tryCatch(
        haseman_elston(
            simulated$genotype, simulated$phenotype, simulated$covariates
        ),
        error = function(e) {
            errors <<- c(errors, paste0("HE: ", conditionMessage(e)))
            data.frame(
                he_h2 = NA_real_, he_se = NA_real_, he_pvalue = NA_real_,
                he_num_snps = NA_integer_, he_converged = FALSE
            )
        }
    )
    p_eff <- tryCatch(
        effective_rank_genotype(simulated$genotype),
        error = function(e) {
            errors <<- c(errors, paste0("p_eff: ", conditionMessage(e)))
            NA_real_
        }
    )

    if (identical(scenario$feature_mode[[1L]], "development_augment")) {
        source_path <- file.path(
            development_bslmm_root,
            sprintf("scenario-%07d.tsv", scenario$scenario_id[[1L]])
        )
        if (!file.exists(source_path)) {
            stop("Missing development BSLMM result: ", source_path)
        }
        bslmm_source <- read_tsv(source_path)
        if (nrow(bslmm_source) != 1L) stop("BSLMM source must have one row")
        assert_source_identity(bslmm_source, scenario, simulated)
        bslmm <- list(
            converged = bslmm_source$bslmm_converged[[1L]] %in% TRUE,
            exit_status = as.integer(bslmm_source$bslmm_exit_status[[1L]]),
            error = bslmm_source$bslmm_error[[1L]],
            elapsed_sec = as.numeric(bslmm_source$bslmm_elapsed_sec[[1L]]),
            n_mcmc = as.integer(bslmm_source$bslmm_n_mcmc[[1L]]),
            pve_mean = as.numeric(bslmm_source$bslmm_pve[[1L]]),
            pve_median = as.numeric(bslmm_source$bslmm_pve_median[[1L]]),
            pve_q025 = as.numeric(bslmm_source$bslmm_pve_q025[[1L]]),
            pve_q975 = as.numeric(bslmm_source$bslmm_pve_q975[[1L]]),
            h_mean = as.numeric(bslmm_source$bslmm_h_mean[[1L]])
        )
    } else if (identical(scenario$feature_mode[[1L]], "full")) {
        residual_y <- residualize_phenotype(
            simulated$phenotype, simulated$covariates
        )
        bslmm_dir <- file.path(
            work_dir, sprintf("scenario-%07d", scenario$scenario_id[[1L]])
        )
        if (dir.exists(bslmm_dir)) unlink(bslmm_dir, recursive = TRUE)
        bslmm <- fit_bslmm_pve(
            genotype = simulated$genotype,
            phenotype = residual_y,
            work_dir = bslmm_dir,
            gemma_bin = scenario$gemma_bin[[1L]],
            bslmm_mode = as.integer(scenario$bslmm_mode[[1L]]),
            burn_in = as.integer(scenario$bslmm_burn_in[[1L]]),
            sampling = as.integer(scenario$bslmm_sampling[[1L]]),
            rpace = as.integer(scenario$bslmm_rpace[[1L]]),
            seed = scenario$seed[[1L]]
        )
        if (!isTRUE(keep_work)) unlink(bslmm_dir, recursive = TRUE)
    } else {
        stop("Unknown feature_mode: ", scenario$feature_mode[[1L]])
    }

    if (!isTRUE(bslmm$converged)) {
        errors <- c(errors, paste0("BSLMM: ", bslmm$error %||% "not converged"))
    }
    required_values <- c(
        bslmm_pve = bslmm$pve_mean,
        he_h2 = he$he_h2[[1L]],
        rho2_oof = en$rho2_oof[[1L]],
        r2_oof = en$r2_oof[[1L]],
        p_eff = p_eff,
        ld_metric = simulated$ld_metric
    )
    feature_complete <- all(is.finite(required_values)) &&
        isTRUE(en$converged[[1L]]) && isTRUE(he$he_converged[[1L]]) &&
        isTRUE(bslmm$converged)

    data.frame(
        scenario_id = scenario$scenario_id,
        split = scenario$split,
        feature_mode = scenario$feature_mode,
        seed = scenario$seed,
        n = scenario$n,
        num_snps = scenario$num_snps,
        ld_rho = scenario$ld_rho,
        architecture = scenario$architecture,
        true_h2 = scenario$true_h2,
        replicate = scenario$replicate,
        realized_h2 = simulated$realized_h2,
        num_causal = length(simulated$causal_index),
        ld_metric = simulated$ld_metric,
        p_eff = p_eff,
        mean_maf = mean(colMeans(simulated$genotype) / 2),
        bslmm_pve = bslmm$pve_mean,
        bslmm_pve_median = bslmm$pve_median,
        bslmm_pve_q025 = bslmm$pve_q025,
        bslmm_pve_q975 = bslmm$pve_q975,
        bslmm_h_mean = bslmm$h_mean,
        bslmm_converged = isTRUE(bslmm$converged),
        bslmm_exit_status = bslmm$exit_status,
        bslmm_elapsed_sec = bslmm$elapsed_sec,
        bslmm_n_mcmc = bslmm$n_mcmc,
        he_h2 = he$he_h2[[1L]],
        he_se = he$he_se[[1L]],
        he_pvalue = he$he_pvalue[[1L]],
        he_converged = isTRUE(he$he_converged[[1L]]),
        rho2_oof = en$rho2_oof[[1L]],
        r2_oof = en$r2_oof[[1L]],
        covariance_ratio_oof = en$covariance_ratio_oof[[1L]],
        score_variance_ratio_oof = en$score_variance_ratio_oof[[1L]],
        en_converged = isTRUE(en$converged[[1L]]),
        feature_complete = feature_complete,
        computational_failure = !feature_complete,
        feature_error = if (length(errors)) paste(errors, collapse = " | ") else NA_character_,
        stringsAsFactors = FALSE
    )
}

joint_signal_spec <- function(settings) {
    list(
        bslmm_pve = list(
            clip = split_numeric(settings$bslmm_clip),
            knots = split_numeric(settings$bslmm_knots)
        ),
        he_h2 = list(
            clip = split_numeric(settings$he_clip),
            knots = split_numeric(settings$he_knots)
        ),
        rho2_oof = list(
            clip = split_numeric(settings$rho2_clip),
            knots = split_numeric(settings$rho2_knots)
        ),
        r2_oof = list(
            clip = split_numeric(settings$r2_clip),
            knots = split_numeric(settings$r2_knots)
        )
    )
}

make_joint_matrix <- function(data, signal_spec, design_scaler = NULL) {
    signal_columns <- list()
    for (feature in names(signal_spec)) {
        x <- as.numeric(data[[feature]])
        if (any(!is.finite(x))) stop("Nonfinite joint feature: ", feature)
        bounds <- signal_spec[[feature]]$clip
        if (length(bounds) != 2L || bounds[[1L]] >= bounds[[2L]]) {
            stop("Invalid clipping bounds for ", feature)
        }
        x <- pmin(bounds[[2L]], pmax(bounds[[1L]], x))
        for (knot in signal_spec[[feature]]$knots) {
            label <- gsub("-", "m", format(knot, scientific = FALSE, trim = TRUE))
            label <- gsub("\\.", "p", label)
            name <- paste0("signal__", feature, "__k", label)
            signal_columns[[name]] <- pmax(x - knot, 0)
        }
    }
    signal_matrix <- do.call(cbind, signal_columns)
    storage.mode(signal_matrix) <- "double"

    design_raw <- cbind(
        design__log_n = log(as.numeric(data$n)),
        design__log_p = log(as.numeric(data$num_snps)),
        design__log_p_eff = log(as.numeric(data$p_eff)),
        design__ld = as.numeric(data$ld_metric)
    )
    if (any(!is.finite(design_raw))) stop("Nonfinite joint design feature")
    if (is.null(design_scaler)) {
        center <- colMeans(design_raw)
        scale_value <- apply(design_raw, 2L, stats::sd)
        scale_value[!is.finite(scale_value) | scale_value <= 0] <- 1
        design_scaler <- list(center = center, scale = scale_value)
    }
    design_matrix <- sweep(design_raw, 2L, design_scaler$center, "-")
    design_matrix <- sweep(design_matrix, 2L, design_scaler$scale, "/")
    list(
        matrix = cbind(signal_matrix, design_matrix),
        design_scaler = design_scaler
    )
}

joint_case_weights <- function(true_h2) {
    ifelse(abs(true_h2 - 0) < 1e-12, 8,
           ifelse(abs(true_h2 - 0.05) < 1e-12, 6,
                  ifelse(abs(true_h2 - 0.10) < 1e-12, 4, 1)))
}

fit_joint_pve_model <- function(fit_data, calibration_data, settings) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("glmnet is required for the locked joint estimator")
    }
    if (!all(fit_data$feature_complete %in% TRUE) ||
        !all(calibration_data$feature_complete %in% TRUE)) {
        stop("Development features are incomplete; estimator cannot be fit")
    }
    spec <- joint_signal_spec(settings)
    built <- make_joint_matrix(fit_data, spec)
    x <- built$matrix
    signal <- startsWith(colnames(x), "signal__")
    lower <- ifelse(signal, 0, -Inf)
    lambda <- as_num(settings$ridge_lambda, "ridge_lambda")
    fit <- glmnet::glmnet(
        x = x,
        y = fit_data$true_h2,
        family = "gaussian",
        alpha = 0,
        lambda = lambda,
        weights = joint_case_weights(fit_data$true_h2),
        intercept = TRUE,
        standardize = TRUE,
        lower.limits = lower,
        thresh = 1e-10,
        maxit = 100000L
    )
    model <- list(
        family = settings$estimator_family,
        gate_version = "final_joint_pve_v1",
        signal_spec = spec,
        design_scaler = built$design_scaler,
        glmnet_fit = fit,
        lambda = lambda,
        output_lower = as_num(settings$output_lower, "output_lower"),
        output_upper = as_num(settings$output_upper, "output_upper"),
        config = settings,
        conformal_q = NA_real_,
        null_cutoff = NA_real_
    )
    cal_prediction <- predict_joint_pve(model, calibration_data)
    alpha <- as_num(settings$conformal_alpha, "conformal_alpha")
    absolute_residual <- abs(
        calibration_data$true_h2 - cal_prediction$pve_cis_joint_calibrated
    )
    m <- length(absolute_residual)
    q_index <- min(m, ceiling((m + 1) * (1 - alpha)))
    model$conformal_q <- sort(absolute_residual)[[q_index]]
    null_prediction <- cal_prediction$pve_cis_joint_calibrated[
        calibration_data$true_h2 == 0
    ]
    m0 <- length(null_prediction)
    null_alpha <- as_num(settings$null_alpha, "null_alpha")
    null_index <- min(m0, ceiling((m0 + 1) * (1 - null_alpha)))
    model$null_cutoff <- sort(null_prediction)[[null_index]]
    model$development_counts <- c(
        fit = nrow(fit_data), calibration = nrow(calibration_data), null = m0
    )
    model
}

predict_joint_pve <- function(model, newdata) {
    ## A calibrator restored with readRDS() does not itself load glmnet's S3
    ## prediction method into a fresh R session. Load the namespace explicitly
    ## before dispatch; this changes no fitted object, coefficient, or estimate.
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("glmnet is required to apply the frozen joint PVE estimator")
    }
    built <- make_joint_matrix(
        newdata, model$signal_spec, design_scaler = model$design_scaler
    )
    unbounded <- drop(stats::predict(
        model$glmnet_fit, newx = built$matrix, s = model$lambda
    ))
    estimate <- pmin(model$output_upper, pmax(model$output_lower, unbounded))
    lower <- if (is.finite(model$conformal_q)) {
        pmax(model$output_lower, estimate - model$conformal_q)
    } else {
        rep(NA_real_, length(estimate))
    }
    upper <- if (is.finite(model$conformal_q)) {
        pmin(model$output_upper, estimate + model$conformal_q)
    } else {
        rep(NA_real_, length(estimate))
    }
    data.frame(
        pve_cis_joint_unbounded = unbounded,
        pve_cis_joint_calibrated = estimate,
        pve_simulation_reference_lower = lower,
        pve_simulation_reference_upper = upper,
        pve_lower_boundary_hit = unbounded <= model$output_lower,
        pve_upper_boundary_hit = unbounded >= model$output_upper,
        positive_signal = if (is.finite(model$null_cutoff)) {
            estimate > model$null_cutoff
        } else {
            NA
        }
    )
}

mean_or_na_joint <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) mean(x) else NA_real_
}

joint_validation_metrics <- function(data, expected_n) {
    available <- data$feature_complete %in% TRUE &
        is.finite(data$pve_cis_joint_calibrated)
    usable <- data[available, , drop = FALSE]
    error <- usable$pve_cis_joint_calibrated - usable$true_h2
    low <- usable$true_h2 <= 0.10 + 1e-12
    null <- usable$true_h2 == 0
    covered <- usable$true_h2 >= usable$pve_simulation_reference_lower &
        usable$true_h2 <= usable$pve_simulation_reference_upper
    by_level <- aggregate(
        cbind(estimate = usable$pve_cis_joint_calibrated,
              error = error,
              squared_error = error^2),
        by = list(true_h2 = usable$true_h2),
        FUN = mean
    )
    by_level$rmse <- sqrt(by_level$squared_error)
    low_means <- by_level$estimate[match(c(0, 0.05, 0.10), by_level$true_h2)]
    min_low_difference <- if (all(is.finite(low_means))) {
        min(diff(low_means))
    } else {
        NA_real_
    }
    metrics <- c(
        absolute_mean_bias = abs(mean_or_na_joint(error)),
        rmse = sqrt(mean_or_na_joint(error^2)),
        calibration_interval_coverage = mean_or_na_joint(covered),
        null_type1_error = mean_or_na_joint(usable$positive_signal[null]),
        spearman_truth_estimate = suppressWarnings(stats::cor(
            usable$true_h2, usable$pve_cis_joint_calibrated,
            method = "spearman", use = "complete.obs"
        )),
        null_mean_estimated_h2 = mean_or_na_joint(
            usable$pve_cis_joint_calibrated[null]
        ),
        low_h2_absolute_mean_bias = abs(mean_or_na_joint(error[low])),
        low_h2_rmse = sqrt(mean_or_na_joint(error[low]^2)),
        max_absolute_low_h2_level_bias = max(
            abs(by_level$error[by_level$true_h2 <= 0.10 + 1e-12]),
            na.rm = TRUE
        ),
        minimum_low_h2_adjacent_mean_difference = min_low_difference,
        max_absolute_h2_level_bias = max(abs(by_level$error), na.rm = TRUE),
        complete_feature_rate = mean(data$feature_complete %in% TRUE),
        computational_failure_rate = mean(data$computational_failure %in% TRUE),
        reconciliation_rate = nrow(data) / expected_n
    )
    list(metrics = metrics, by_level = by_level)
}

evaluate_joint_criteria <- function(metrics, criteria) {
    value <- unname(metrics[criteria$metric])
    if (any(!is.finite(value))) {
        missing <- criteria$metric[!is.finite(value)]
        stop("Acceptance metrics missing/nonfinite: ", paste(missing, collapse = ", "))
    }
    threshold <- as.numeric(criteria$threshold)
    passed <- mapply(function(x, op, limit) {
        switch(
            op,
            less_than_or_equal = x <= limit,
            greater_than_or_equal = x >= limit,
            greater_than = x > limit,
            less_than = x < limit,
            stop("Unknown acceptance comparison: ", op)
        )
    }, value, criteria$comparison, threshold)
    data.frame(criteria, value = value, passed = as.logical(passed),
               stringsAsFactors = FALSE)
}
