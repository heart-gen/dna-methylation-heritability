## Shared functions for simulation-calibrated elastic-net estimates of local
## SNP-explained methylation variance.
##
## All phenotype-informed preprocessing is performed inside an outer training
## fold. The held-out predictions are the only predictions used for raw
## estimator calculation.

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L || is.na(x) || identical(x, "")) y else x
}

parse_cli <- function(defaults = list()) {
    args <- commandArgs(trailingOnly = TRUE)
    out <- defaults
    i <- 1L
    while (i <= length(args)) {
        token <- args[[i]]
        if (!startsWith(token, "--")) {
            stop("Unexpected positional argument: ", token)
        }
        token <- substring(token, 3L)
        if (grepl("=", token, fixed = TRUE)) {
            pieces <- strsplit(token, "=", fixed = TRUE)[[1L]]
            key <- pieces[[1L]]
            value <- paste(pieces[-1L], collapse = "=")
        } else {
            key <- token
            if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
                value <- "TRUE"
            } else {
                i <- i + 1L
                value <- args[[i]]
            }
        }
        key <- gsub("-", "_", key, fixed = TRUE)
        out[[key]] <- value
        i <- i + 1L
    }
    out
}

as_int <- function(x, name) {
    value <- suppressWarnings(as.integer(x))
    if (length(value) != 1L || is.na(value)) stop(name, " must be an integer")
    value
}

as_num <- function(x, name) {
    value <- suppressWarnings(as.numeric(x))
    if (length(value) != 1L || is.na(value)) stop(name, " must be numeric")
    value
}

as_bool <- function(x, name) {
    value <- tolower(as.character(x))
    if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
    if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
    stop(name, " must be TRUE or FALSE")
}

split_numeric <- function(x) {
    values <- suppressWarnings(as.numeric(strsplit(as.character(x), ",", fixed = TRUE)[[1L]]))
    if (anyNA(values)) stop("Could not parse numeric vector: ", x)
    values
}

read_tsv <- function(path, ...) {
    if (!file.exists(path)) stop("Input file does not exist: ", path)
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

write_tsv <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
    write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    if (!file.rename(tmp, path)) stop("Could not atomically write: ", path)
    invisible(path)
}

make_balanced_folds <- function(n, k, seed) {
    if (n < 6L) stop("At least six samples are required")
    k <- max(2L, min(as.integer(k), n))
    set.seed(seed)
    assignment <- rep(seq_len(k), length.out = n)
    sample(assignment, n, replace = FALSE)
}

impute_from_training <- function(x_train, x_test) {
    means <- colMeans(x_train, na.rm = TRUE)
    means[!is.finite(means)] <- 0
    for (j in seq_len(ncol(x_train))) {
        miss_train <- is.na(x_train[, j])
        miss_test <- is.na(x_test[, j])
        if (any(miss_train)) x_train[miss_train, j] <- means[[j]]
        if (any(miss_test)) x_test[miss_test, j] <- means[[j]]
    }
    list(train = x_train, test = x_test)
}

adjust_phenotype_in_fold <- function(y_train, y_test, covar_train = NULL,
                                     covar_test = NULL) {
    if (is.null(covar_train) || ncol(covar_train) == 0L) {
        train_adjusted <- y_train - mean(y_train)
        test_adjusted <- y_test - mean(y_train)
    } else {
        x_train <- cbind(`(Intercept)` = 1, as.matrix(covar_train))
        x_test <- cbind(`(Intercept)` = 1, as.matrix(covar_test))
        fit <- lm.fit(x = x_train, y = y_train)
        coefficients <- fit$coefficients
        coefficients[!is.finite(coefficients)] <- 0
        train_adjusted <- y_train - drop(x_train %*% coefficients)
        test_adjusted <- y_test - drop(x_test %*% coefficients)
    }
    center <- mean(train_adjusted)
    scale_value <- stats::sd(train_adjusted)
    if (!is.finite(scale_value) || scale_value <= 1e-10) {
        stop("Training phenotype has no residual variance")
    }
    list(
        train = (train_adjusted - center) / scale_value,
        test = (test_adjusted - center) / scale_value,
        center = center,
        scale = scale_value
    )
}

screen_training_snps <- function(x_train, y_train, max_features) {
    p <- ncol(x_train)
    if (p <= max_features) return(seq_len(p))
    x_centered <- sweep(x_train, 2L, colMeans(x_train), FUN = "-")
    denom <- sqrt(colSums(x_centered^2) * sum((y_train - mean(y_train))^2))
    score <- abs(drop(crossprod(x_centered, y_train - mean(y_train))) / denom)
    score[!is.finite(score)] <- -Inf
    order(score, decreasing = TRUE)[seq_len(max_features)]
}

fit_inner_elastic_net <- function(x, y, alpha_grid, inner_folds, seed,
                                  lambda_rule = "lambda.1se") {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("The glmnet package is required")
    }
    inner_folds <- max(3L, min(as.integer(inner_folds), nrow(x)))
    foldid <- make_balanced_folds(nrow(x), inner_folds, seed)
    candidates <- vector("list", length(alpha_grid))
    scores <- rep(Inf, length(alpha_grid))

    for (i in seq_along(alpha_grid)) {
        alpha <- alpha_grid[[i]]
        fit <- tryCatch(
            glmnet::cv.glmnet(
                x = x,
                y = y,
                alpha = alpha,
                foldid = foldid,
                family = "gaussian",
                standardize = TRUE,
                intercept = TRUE,
                type.measure = "mse",
                keep = FALSE,
                parallel = FALSE
            ),
            error = function(e) NULL
        )
        candidates[[i]] <- fit
        if (!is.null(fit)) {
            lambda_value <- fit[[lambda_rule]]
            lambda_index <- which.min(abs(fit$lambda - lambda_value))
            scores[[i]] <- fit$cvm[[lambda_index]]
        }
    }
    best <- which.min(scores)
    if (!length(best) || !is.finite(scores[[best]]) || is.null(candidates[[best]])) {
        stop("All inner elastic-net fits failed")
    }
    list(
        fit = candidates[[best]],
        alpha = alpha_grid[[best]],
        lambda = candidates[[best]][[lambda_rule]],
        inner_mse = scores[[best]],
        lambda_rule = lambda_rule
    )
}

safe_ratio <- function(numerator, denominator) {
    if (!is.finite(denominator) || abs(denominator) <= 1e-12) return(NA_real_)
    numerator / denominator
}

squared_prediction_correlation <- function(observed, predicted,
                                           variance_tolerance = 1e-12) {
    prediction_variance <- stats::var(predicted)
    ## Test the null-model condition before computing correlation. A vector of
    ## floating-point noise can have a finite correlation despite carrying no
    ## usable predictive variation.
    if (is.finite(prediction_variance) &&
        prediction_variance <= variance_tolerance) {
        return(0)
    }
    rho <- suppressWarnings(stats::cor(observed, predicted))
    if (is.finite(rho)) rho^2 else NA_real_
}

crossfit_elastic_net <- function(genotype, phenotype, covariates = NULL,
                                 outer_folds = 5L, outer_repeats = 5L,
                                 inner_folds = 5L,
                                 alpha_grid = c(0.05, 0.25, 0.5, 0.75, 1),
                                 lambda_rule = "lambda.1se",
                                 max_features = 2000L,
                                 seed = 20250805L,
                                 keep_predictions = TRUE) {
    genotype <- as.matrix(genotype)
    storage.mode(genotype) <- "double"
    phenotype <- as.numeric(phenotype)
    if (nrow(genotype) != length(phenotype)) stop("Genotype and phenotype sample counts differ")
    if (any(!is.finite(phenotype))) stop("Phenotype contains missing or non-finite values")
    if (!is.null(covariates)) {
        covariates <- as.matrix(covariates)
        storage.mode(covariates) <- "double"
        if (nrow(covariates) != length(phenotype)) stop("Covariate and phenotype sample counts differ")
        if (any(!is.finite(covariates))) stop("Covariates contain missing or non-finite values")
    }
    n <- nrow(genotype)
    p <- ncol(genotype)
    if (p < 1L) stop("No SNPs supplied")
    outer_folds <- max(2L, min(as.integer(outer_folds), n))
    outer_repeats <- as.integer(outer_repeats)
    prediction_sum <- numeric(n)
    adjusted_sum <- numeric(n)
    prediction_count <- integer(n)
    fold_records <- list()
    record_index <- 1L

    for (repeat_id in seq_len(outer_repeats)) {
        foldid <- make_balanced_folds(n, outer_folds, seed + repeat_id * 1009L)
        for (fold in seq_len(outer_folds)) {
            test_index <- which(foldid == fold)
            train_index <- which(foldid != fold)
            x_train <- genotype[train_index, , drop = FALSE]
            x_test <- genotype[test_index, , drop = FALSE]
            imputed <- impute_from_training(x_train, x_test)
            x_train <- imputed$train
            x_test <- imputed$test

            train_variance <- apply(x_train, 2L, stats::var)
            keep <- which(is.finite(train_variance) & train_variance > 1e-8)
            if (!length(keep)) stop("No polymorphic SNPs in outer training fold")
            x_train <- x_train[, keep, drop = FALSE]
            x_test <- x_test[, keep, drop = FALSE]

            covar_train <- if (is.null(covariates)) NULL else covariates[train_index, , drop = FALSE]
            covar_test <- if (is.null(covariates)) NULL else covariates[test_index, , drop = FALSE]
            adjusted <- adjust_phenotype_in_fold(
                phenotype[train_index], phenotype[test_index], covar_train, covar_test
            )

            selected <- screen_training_snps(x_train, adjusted$train, max_features)
            x_train_selected <- x_train[, selected, drop = FALSE]
            x_test_selected <- x_test[, selected, drop = FALSE]
            tuned <- fit_inner_elastic_net(
                x_train_selected,
                adjusted$train,
                alpha_grid = alpha_grid,
                inner_folds = inner_folds,
                seed = seed + repeat_id * 1009L + fold * 9173L,
                lambda_rule = lambda_rule
            )
            prediction <- drop(stats::predict(
                tuned$fit,
                newx = x_test_selected,
                s = tuned$lambda
            ))
            coefficients <- as.numeric(stats::coef(tuned$fit, s = tuned$lambda))[-1L]

            prediction_sum[test_index] <- prediction_sum[test_index] + prediction
            adjusted_sum[test_index] <- adjusted_sum[test_index] + adjusted$test
            prediction_count[test_index] <- prediction_count[test_index] + 1L
            fold_records[[record_index]] <- data.frame(
                repeat_id = repeat_id,
                fold = fold,
                n_train = length(train_index),
                n_test = length(test_index),
                snps_polymorphic = length(keep),
                snps_screened = length(selected),
                snps_nonzero = sum(coefficients != 0),
                alpha = tuned$alpha,
                lambda = tuned$lambda,
                inner_mse = tuned$inner_mse,
                fold_score_variance_ratio = safe_ratio(
                    stats::var(prediction), stats::var(adjusted$test)
                )
            )
            record_index <- record_index + 1L
        }
    }
    if (any(prediction_count != outer_repeats)) {
        stop("Every sample must receive exactly one held-out prediction per repeat")
    }
    prediction <- prediction_sum / prediction_count
    adjusted_phenotype <- adjusted_sum / prediction_count
    phenotype_variance <- stats::var(adjusted_phenotype)
    centered_sst <- sum((adjusted_phenotype - mean(adjusted_phenotype))^2)
    r2_oof <- 1 - safe_ratio(sum((adjusted_phenotype - prediction)^2), centered_sst)
    prediction_variance <- stats::var(prediction)
    rho2_oof <- squared_prediction_correlation(adjusted_phenotype, prediction)
    covariance_ratio <- safe_ratio(stats::cov(adjusted_phenotype, prediction), phenotype_variance)
    score_variance_ratio <- safe_ratio(prediction_variance, phenotype_variance)
    calibration_slope <- safe_ratio(stats::cov(adjusted_phenotype, prediction), prediction_variance)
    fold_details <- do.call(rbind, fold_records)

    result <- list(
        metrics = data.frame(
            n = n,
            num_snps = p,
            outer_folds = outer_folds,
            outer_repeats = outer_repeats,
            inner_folds = inner_folds,
            max_features = max_features,
            r2_oof = r2_oof,
            rho2_oof = rho2_oof,
            covariance_ratio_oof = covariance_ratio,
            score_variance_ratio_oof = score_variance_ratio,
            calibration_slope_oof = calibration_slope,
            mean_fold_score_variance_ratio = mean(fold_details$fold_score_variance_ratio, na.rm = TRUE),
            mean_nonzero_snps = mean(fold_details$snps_nonzero),
            converged = TRUE
        ),
        folds = fold_details
    )
    if (keep_predictions) {
        result$predictions <- data.frame(
            sample_index = seq_len(n),
            adjusted_phenotype = adjusted_phenotype,
            oof_prediction = prediction,
            prediction_repeats = prediction_count
        )
    }
    result
}

haseman_elston <- function(genotype, phenotype, covariates = NULL) {
    ## Method-of-moments sensitivity estimator. This avoids REML optimization,
    ## but it can be negative or exceed one and is expected to be imprecise at
    ## the study sample sizes. Values are therefore never clipped.
    genotype <- as.matrix(genotype)
    storage.mode(genotype) <- "double"
    phenotype <- as.numeric(phenotype)
    if (nrow(genotype) != length(phenotype)) stop("Genotype and phenotype sample counts differ")
    if (anyNA(genotype)) {
        means <- colMeans(genotype, na.rm = TRUE)
        means[!is.finite(means)] <- 0
        for (j in seq_len(ncol(genotype))) {
            missing <- is.na(genotype[, j])
            if (any(missing)) genotype[missing, j] <- means[[j]]
        }
    }
    variances <- apply(genotype, 2L, stats::var)
    keep <- is.finite(variances) & variances > 1e-8
    genotype <- genotype[, keep, drop = FALSE]
    if (ncol(genotype) < 2L) {
        return(data.frame(
            he_h2 = NA_real_, he_se = NA_real_, he_pvalue = NA_real_,
            he_num_snps = ncol(genotype), he_converged = FALSE
        ))
    }
    standardized <- scale(genotype)
    grm <- tcrossprod(standardized) / ncol(standardized)
    if (is.null(covariates) || ncol(as.matrix(covariates)) == 0L) {
        residual <- phenotype - mean(phenotype)
    } else {
        design <- cbind(`(Intercept)` = 1, as.matrix(covariates))
        coefficients <- lm.fit(design, phenotype)$coefficients
        coefficients[!is.finite(coefficients)] <- 0
        residual <- phenotype - drop(design %*% coefficients)
    }
    residual <- as.numeric(scale(residual))
    upper <- upper.tri(grm)
    response <- outer(residual, residual)[upper]
    relatedness <- grm[upper]
    design <- cbind(`(Intercept)` = 1, relatedness = relatedness)
    fit <- lm.fit(design, response)
    degrees_freedom <- length(response) - fit$rank
    coefficient <- fit$coefficients[[2L]]
    if (!is.finite(coefficient) || degrees_freedom <= 0L) {
        return(data.frame(
            he_h2 = NA_real_, he_se = NA_real_, he_pvalue = NA_real_,
            he_num_snps = ncol(genotype), he_converged = FALSE
        ))
    }
    sigma2 <- sum(fit$residuals^2) / degrees_freedom
    covariance <- tryCatch(
        sigma2 * chol2inv(qr.R(fit$qr)[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE]),
        error = function(e) matrix(NA_real_, nrow = 2L, ncol = 2L)
    )
    standard_error <- sqrt(covariance[2L, 2L])
    statistic <- safe_ratio(coefficient, standard_error)
    pvalue <- if (is.finite(statistic)) {
        2 * stats::pt(abs(statistic), df = degrees_freedom, lower.tail = FALSE)
    } else {
        NA_real_
    }
    data.frame(
        he_h2 = coefficient,
        he_se = standard_error,
        he_pvalue = pvalue,
        he_num_snps = ncol(genotype),
        he_converged = is.finite(coefficient) && is.finite(standard_error)
    )
}

simulate_ar1_genotypes <- function(n, p, ld_rho, maf_min = 0.05,
                                   maf_max = 0.5) {
    if (abs(ld_rho) >= 1) stop("ld_rho must be between -1 and 1")
    latent <- matrix(0, nrow = n, ncol = p)
    latent[, 1L] <- stats::rnorm(n)
    innovation_scale <- sqrt(1 - ld_rho^2)
    if (p > 1L) {
        for (j in 2:p) {
            latent[, j] <- ld_rho * latent[, j - 1L] +
                innovation_scale * stats::rnorm(n)
        }
    }
    maf <- stats::runif(p, maf_min, maf_max)
    genotype <- matrix(0, nrow = n, ncol = p)
    for (j in seq_len(p)) {
        q <- 1 - maf[[j]]
        threshold_0 <- stats::qnorm(q^2)
        threshold_1 <- stats::qnorm(q^2 + 2 * maf[[j]] * q)
        genotype[, j] <- ifelse(
            latent[, j] <= threshold_0, 0,
            ifelse(latent[, j] <= threshold_1, 1, 2)
        )
    }
    colnames(genotype) <- paste0("snp_", seq_len(p))
    genotype
}

adjacent_ld_metric <- function(genotype, max_pairs = 200L) {
    p <- ncol(genotype)
    if (p < 2L) return(0)
    indices <- unique(round(seq(1, p - 1L, length.out = min(max_pairs, p - 1L))))
    values <- vapply(indices, function(j) {
        value <- suppressWarnings(stats::cor(genotype[, j], genotype[, j + 1L], use = "pairwise.complete.obs"))
        if (is.finite(value)) value^2 else NA_real_
    }, numeric(1L))
    stats::median(values, na.rm = TRUE)
}

simulate_locus <- function(n, p, ld_rho, h2, architecture,
                           causal_count = NA_integer_, causal_fraction = NA_real_) {
    genotype <- simulate_ar1_genotypes(n, p, ld_rho)
    if (architecture == "sparse") {
        number_causal <- if (is.na(causal_count)) min(5L, p) else min(causal_count, p)
    } else if (architecture == "oligogenic") {
        fraction <- if (is.na(causal_fraction)) 0.01 else causal_fraction
        number_causal <- min(p, max(10L, ceiling(p * fraction)))
    } else if (architecture == "polygenic") {
        fraction <- if (is.na(causal_fraction)) 1 else causal_fraction
        number_causal <- min(p, max(1L, ceiling(p * fraction)))
    } else {
        stop("Unknown architecture: ", architecture)
    }
    causal_index <- if (h2 > 0) sample(seq_len(p), number_causal) else integer()
    beta <- numeric(p)
    if (length(causal_index)) beta[causal_index] <- stats::rnorm(length(causal_index))
    genetic_value <- drop(genotype %*% beta)
    if (h2 > 0 && stats::var(genetic_value) > 0) {
        genetic_value <- as.numeric(scale(genetic_value)) * sqrt(h2)
    } else {
        genetic_value[] <- 0
    }
    noise <- stats::rnorm(n)
    if (stats::var(genetic_value) > 0) {
        noise <- stats::residuals(stats::lm(noise ~ genetic_value))
    }
    noise <- as.numeric(scale(noise)) * sqrt(1 - h2)
    residual_phenotype <- genetic_value + noise

    age <- as.numeric(scale(stats::rnorm(n, 50, 10)))
    sex <- stats::rbinom(n, 1, 0.5)
    batch <- factor(sample(seq_len(3L), n, replace = TRUE))
    covariates <- stats::model.matrix(~ age + sex + batch)[, -1L, drop = FALSE]
    covariate_beta <- seq(0.15, 0.05, length.out = ncol(covariates))
    phenotype <- residual_phenotype + drop(covariates %*% covariate_beta)

    list(
        genotype = genotype,
        phenotype = phenotype,
        covariates = covariates,
        genetic_value = genetic_value,
        residual_phenotype = residual_phenotype,
        beta = beta,
        causal_index = causal_index,
        realized_h2 = safe_ratio(stats::var(genetic_value), stats::var(residual_phenotype)),
        ld_metric = adjacent_ld_metric(genotype)
    )
}

fit_isotonic_curve <- function(raw, truth) {
    keep <- is.finite(raw) & is.finite(truth)
    raw <- raw[keep]
    truth <- truth[keep]
    if (length(raw) < 10L || length(unique(raw)) < 3L) {
        stop("Insufficient observations to fit isotonic calibration curve")
    }
    order_index <- order(raw, truth)
    fit <- stats::isoreg(raw[order_index], truth[order_index])
    curve <- aggregate(fit$yf, by = list(raw = fit$x), FUN = mean)
    names(curve)[[2L]] <- "calibrated"
    curve$calibrated <- pmin(1, pmax(0, curve$calibrated))
    fitted <- stats::approx(curve$raw, curve$calibrated, xout = raw,
                            rule = 2, ties = mean)$y
    residual <- truth - fitted
    list(
        curve = curve,
        residual_quantiles = stats::quantile(
            residual, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE
        ),
        raw_range = range(raw),
        truth_range = range(truth),
        n = length(raw)
    )
}

make_forward_features <- function(data, boundary_aware = FALSE) {
    required <- c(
        "rho2_oof", "r2_oof", "covariance_ratio_oof",
        "score_variance_ratio_oof", "mean_nonzero_snps"
    )
    missing <- setdiff(required, names(data))
    if (length(missing)) {
        stop("Forward calibration input is missing: ", paste(missing, collapse = ", "))
    }
    features <- data.frame(
        sqrt_rho = sqrt(pmax(as.numeric(data$rho2_oof), 0)),
        r2 = as.numeric(data$r2_oof),
        covariance = as.numeric(data$covariance_ratio_oof),
        sqrt_score_variance = sqrt(pmax(as.numeric(data$score_variance_ratio_oof), 0)),
        log_nonzero_snps = log1p(pmax(as.numeric(data$mean_nonzero_snps), 0)),
        stringsAsFactors = FALSE
    )
    if (isTRUE(boundary_aware)) {
        he <- if ("he_h2" %in% names(data)) as.numeric(data$he_h2) else rep(NA_real_, nrow(data))
        features$near_ceiling <- pmax(0, features$sqrt_rho - 0.5)^2
        features$he_excess <- pmax(0, he - 0.6)
        features$high_signal_flag <- as.numeric(features$sqrt_rho >= 0.7)
    }
    features
}

forward_design_id <- function(n, p, ld_rho) {
    paste0("n", n, "_p", p, "_ld", format(ld_rho, trim = TRUE))
}

fit_forward_scalers <- function(data, boundary_aware = FALSE) {
    required <- c("n", "num_snps", "ld_rho", "ld_metric")
    missing <- setdiff(required, names(data))
    if (length(missing)) stop("Calibration design is missing: ", paste(missing, collapse = ", "))
    features <- make_forward_features(data, boundary_aware = boundary_aware)
    ids <- forward_design_id(data$n, data$num_snps, data$ld_rho)
    split_index <- split(seq_len(nrow(data)), ids)
    scalers <- lapply(names(split_index), function(id) {
        index <- split_index[[id]]
        center <- vapply(features[index, , drop = FALSE], function(x) {
            mean(x[is.finite(x)])
        }, numeric(1L))
        scale_value <- vapply(features[index, , drop = FALSE], function(x) {
            stats::sd(x[is.finite(x)])
        }, numeric(1L))
        center[!is.finite(center)] <- 0
        scale_value[!is.finite(scale_value) | scale_value <= 1e-10] <- 1
        list(
            id = id,
            n_center = stats::median(data$n[index]),
            p_center = stats::median(data$num_snps[index]),
            ld_center = stats::median(data$ld_metric[index], na.rm = TRUE),
            ld_rho_design = unique(data$ld_rho[index])[[1L]],
            feature_center = center,
            feature_scale = scale_value,
            boundary_aware = isTRUE(boundary_aware)
        )
    })
    names(scalers) <- vapply(scalers, `[[`, character(1L), "id")
    scalers
}

transform_forward_features <- function(data, scalers) {
    required <- c("n", "num_snps", "ld_metric")
    missing <- setdiff(required, names(data))
    if (length(missing)) stop("Forward prediction design is missing: ", paste(missing, collapse = ", "))
    boundary_aware <- isTRUE(scalers[[1L]]$boundary_aware)
    features <- make_forward_features(data, boundary_aware = boundary_aware)
    selected <- character(nrow(data))
    distance <- numeric(nrow(data))
    for (i in seq_len(nrow(data))) {
        exact <- character()
        if ("ld_rho" %in% names(data) && is.finite(data$ld_rho[[i]])) {
            candidate <- forward_design_id(data$n[[i]], data$num_snps[[i]], data$ld_rho[[i]])
            if (candidate %in% names(scalers)) exact <- candidate
        }
        if (length(exact)) {
            id <- exact[[1L]]
            distances <- calibration_distance(
                data$n[[i]], data$num_snps[[i]], data$ld_metric[[i]], scalers[[id]]
            )
        } else {
            all_distances <- vapply(scalers, function(x) {
                calibration_distance(data$n[[i]], data$num_snps[[i]], data$ld_metric[[i]], x)
            }, numeric(1L))
            id <- names(which.min(all_distances))[[1L]]
            distances <- all_distances[[id]]
        }
        selected[[i]] <- id
        distance[[i]] <- distances
        scaler <- scalers[[id]]
        for (feature in names(features)) {
            value <- features[[feature]][[i]]
            if (!is.finite(value)) value <- scaler$feature_center[[feature]]
            features[[feature]][[i]] <-
                (value - scaler$feature_center[[feature]]) / scaler$feature_scale[[feature]]
        }
    }
    design <- lapply(selected, function(id) scalers[[id]])
    features$n_factor <- factor(
        vapply(design, function(x) as.character(x$n_center), character(1L)),
        levels = sort(unique(vapply(scalers, function(x) as.character(x$n_center), character(1L))))
    )
    features$p_factor <- factor(
        vapply(design, function(x) as.character(x$p_center), character(1L)),
        levels = sort(unique(vapply(scalers, function(x) as.character(x$p_center), character(1L))))
    )
    features$ld_factor <- factor(
        vapply(design, function(x) as.character(x$ld_rho_design), character(1L)),
        levels = sort(unique(vapply(scalers, function(x) as.character(x$ld_rho_design), character(1L))))
    )
    attr(features, "calibration_stratum") <- selected
    attr(features, "calibration_distance") <- distance
    features
}

fit_forward_regression <- function(data, boundary_aware = FALSE) {
    if (!"true_h2" %in% names(data)) stop("true_h2 is required to fit forward calibration")
    scalers <- fit_forward_scalers(data, boundary_aware = boundary_aware)
    transformed <- transform_forward_features(data, scalers)
    transformed$true_h2 <- data$true_h2
    core_terms <- "sqrt_rho + r2 + covariance + sqrt_score_variance + log_nonzero_snps"
    if (isTRUE(boundary_aware)) {
        core_terms <- paste(
            core_terms, "+ near_ceiling + he_excess + high_signal_flag"
        )
    }
    if (nrow(data) >= 100L &&
        length(unique(data$n)) > 1L &&
        length(unique(data$num_snps)) > 1L &&
        length(unique(data$ld_rho)) > 1L) {
        formula <- stats::as.formula(paste(
            "true_h2 ~ (", core_terms, ") * (n_factor + p_factor + ld_factor)"
        ))
    } else {
        ## Small deterministic smoke grids cannot support all design interactions.
        formula <- stats::as.formula(paste("true_h2 ~", core_terms))
    }
    fit <- stats::lm(formula, data = transformed)
    list(
        fit = fit,
        scalers = scalers,
        formula = paste(deparse(formula), collapse = " "),
        boundary_aware = isTRUE(boundary_aware)
    )
}

predict_forward_regression <- function(model, newdata) {
    transformed <- transform_forward_features(newdata, model$scalers)
    estimate <- suppressWarnings(as.numeric(stats::predict(model$fit, newdata = transformed)))
    list(
        estimate = estimate,
        calibration_stratum = attr(transformed, "calibration_stratum"),
        calibration_distance = attr(transformed, "calibration_distance")
    )
}

fit_affine_level_debiasing <- function(score, truth) {
    keep <- is.finite(score) & is.finite(truth)
    level_means <- aggregate(
        score[keep], by = list(true_h2 = truth[keep]), FUN = mean
    )
    names(level_means)[[2L]] <- "forward_score"
    if (nrow(level_means) < 2L || length(unique(level_means$forward_score)) < 2L) {
        stop("At least two distinguishable h2 levels are required for affine debiasing")
    }
    list(
        fit = stats::lm(true_h2 ~ forward_score, data = level_means),
        level_means = level_means
    )
}

predict_affine_level_debiasing <- function(model, score) {
    as.numeric(stats::predict(model$fit, newdata = data.frame(forward_score = score)))
}

fit_two_part_rho_split <- function(data, threshold) {
    if (!is.finite(threshold) || threshold <= 0 || threshold >= 1) {
        stop("two-part rho split threshold must be in (0, 1)")
    }
    features <- make_forward_features(data)
    low_index <- features$sqrt_rho <= threshold
    high_index <- !low_index
    if (sum(low_index) < 20L || sum(high_index) < 20L) {
        stop("Insufficient observations on each side of rho split threshold ", threshold)
    }
    list(
        family = "two_part_rho_split",
        threshold = threshold,
        low_model = fit_forward_regression(data[low_index, , drop = FALSE]),
        high_model = fit_forward_regression(data[high_index, , drop = FALSE])
    )
}

predict_two_part_rho_split <- function(model, newdata) {
    features <- make_forward_features(newdata)
    estimate <- numeric(nrow(newdata))
    stratum <- character(nrow(newdata))
    distance <- numeric(nrow(newdata))
    low_index <- features$sqrt_rho <= model$threshold
    high_index <- !low_index
    if (any(low_index)) {
        low <- predict_forward_regression(model$low_model, newdata[low_index, , drop = FALSE])
        estimate[low_index] <- low$estimate
        stratum[low_index] <- low$calibration_stratum
        distance[low_index] <- low$calibration_distance
    }
    if (any(high_index)) {
        high <- predict_forward_regression(model$high_model, newdata[high_index, , drop = FALSE])
        estimate[high_index] <- high$estimate
        stratum[high_index] <- high$calibration_stratum
        distance[high_index] <- high$calibration_distance
    }
    list(
        estimate = estimate,
        calibration_stratum = stratum,
        calibration_distance = distance
    )
}

fit_two_part_hurdle_high <- function(data, high_threshold = 0.8) {
    if (!is.finite(high_threshold) || high_threshold <= 0 || high_threshold >= 1) {
        stop("hurdle high threshold must be in (0, 1)")
    }
    high_index <- data$true_h2 >= high_threshold
    interior_index <- !high_index
    if (sum(high_index) < 20L || sum(interior_index) < 20L) {
        stop("Insufficient observations for hurdle threshold ", high_threshold)
    }
    gate_scalers <- fit_forward_scalers(data, boundary_aware = TRUE)
    gate_features <- transform_forward_features(data, gate_scalers)
    gate_features$high <- as.integer(high_index)
    gate <- stats::glm(
        high ~ sqrt_rho + r2 + covariance + sqrt_score_variance +
            log_nonzero_snps + near_ceiling + he_excess + high_signal_flag +
            n_factor + p_factor + ld_factor,
        data = gate_features,
        family = stats::binomial()
    )
    list(
        family = "two_part_hurdle_high",
        high_threshold = high_threshold,
        gate = gate,
        gate_scalers = gate_scalers,
        interior_model = fit_forward_regression(data[interior_index, , drop = FALSE]),
        high_model = fit_forward_regression(data[high_index, , drop = FALSE])
    )
}

predict_two_part_hurdle_high <- function(model, newdata) {
    gate_features <- transform_forward_features(newdata, model$gate_scalers)
    p_high <- as.numeric(stats::predict(model$gate, newdata = gate_features, type = "response"))
    interior <- predict_forward_regression(model$interior_model, newdata)
    high <- predict_forward_regression(model$high_model, newdata)
    estimate <- (1 - p_high) * interior$estimate + p_high * high$estimate
    list(
        estimate = estimate,
        calibration_stratum = interior$calibration_stratum,
        calibration_distance = interior$calibration_distance,
        p_high = p_high
    )
}

resolve_calibration_upper_bound <- function(rule, max_simulated_h2) {
    if (identical(rule, "one")) return(1)
    if (identical(rule, "max_simulated_h2")) return(max_simulated_h2)
    stop("Unknown upper_bound_rule: ", rule)
}

fit_candidate_forward_score <- function(candidate_id, family, parameter, data) {
    if (identical(family, "forward_hybrid")) {
        return(list(
            candidate_id = candidate_id,
            family = family,
            forward_model = fit_forward_regression(data)
        ))
    }
    if (identical(family, "boundary_feature_hybrid")) {
        return(list(
            candidate_id = candidate_id,
            family = family,
            forward_model = fit_forward_regression(data, boundary_aware = TRUE)
        ))
    }
    if (identical(family, "two_part_rho_split")) {
        return(c(
            list(candidate_id = candidate_id),
            fit_two_part_rho_split(data, as.numeric(parameter))
        ))
    }
    if (identical(family, "two_part_hurdle_high")) {
        return(c(
            list(candidate_id = candidate_id),
            fit_two_part_hurdle_high(data, as.numeric(parameter))
        ))
    }
    stop("Unknown candidate family: ", family)
}

predict_candidate_forward_score <- function(candidate_model, newdata) {
    family <- candidate_model$family
    if (family %in% c("forward_hybrid", "boundary_feature_hybrid")) {
        return(predict_forward_regression(candidate_model$forward_model, newdata))
    }
    if (identical(family, "two_part_rho_split")) {
        return(predict_two_part_rho_split(candidate_model, newdata))
    }
    if (identical(family, "two_part_hurdle_high")) {
        return(predict_two_part_hurdle_high(candidate_model, newdata))
    }
    stop("Unknown candidate family: ", family)
}

score_hybrid_on_tune <- function(forward_h2, he_h2, truth, upper_bound, weight,
                                 low_h2_max = 0.10) {
    estimate <- pmin(upper_bound, (1 - weight) * forward_h2 + weight * he_h2)
    error <- estimate - truth
    level <- aggregate(error, by = list(true_h2 = truth), FUN = mean)
    low <- truth <= low_h2_max
    low_level <- level[level$true_h2 <= low_h2_max, , drop = FALSE]
    if (!any(low) || !nrow(low_level)) {
        stop("No simulations fall in the configured low-h2 range")
    }
    null <- estimate[truth == 0]
    null_se <- if (length(null) > 1L) {
        stats::sd(null) / sqrt(length(null))
    } else {
        NA_real_
    }
    data.frame(
        he_weight = weight,
        null_mean_estimated_h2 = mean(null),
        null_mean_standard_error = null_se,
        null_mean_upper_95 = mean(null) + stats::qnorm(0.95) * null_se,
        low_h2_absolute_mean_bias = abs(mean(error[low])),
        max_absolute_low_h2_level_bias = max(abs(low_level$x)),
        max_absolute_h2_level_bias = max(abs(level$x)),
        rmse = sqrt(mean(error^2)),
        absolute_mean_bias = abs(mean(error)),
        spearman_truth_estimate = suppressWarnings(stats::cor(
            estimate, truth, method = "spearman"
        )),
        stringsAsFactors = FALSE
    )
}

finite_sample_upper_threshold <- function(values, alpha = 0.05) {
    ## Split-conformal upper cutoff. For an exchangeable future null score,
    ## P(score_new > cutoff) <= alpha in finite samples when the requested
    ## order statistic exists. This is deliberately conservative for small n.
    values <- sort(values[is.finite(values)])
    if (!length(values)) stop("No finite null values supplied")
    if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
        stop("alpha must be strictly between zero and one")
    }
    order_index <- ceiling((length(values) + 1) * (1 - alpha))
    if (order_index > length(values)) {
        return(list(
            threshold = Inf,
            order_index = order_index,
            n = length(values),
            attainable_alpha = 0
        ))
    }
    list(
        threshold = values[[order_index]],
        order_index = order_index,
        n = length(values),
        attainable_alpha = (length(values) + 1 - order_index) /
            (length(values) + 1)
    )
}

calibration_distance <- function(n, p, ld_metric, stratum) {
    abs(log(n / stratum$n_center)) +
        abs(log(p / stratum$p_center)) +
        3 * abs(ld_metric - stratum$ld_center)
}

predict_calibration <- function(calibration_model, newdata,
                                raw_metric = calibration_model$raw_metric) {
    required <- c("n", "num_snps", "ld_metric", raw_metric)
    hybrid_versions <- c("forward_hybrid_v1", "candidate_hybrid_v1")
    if (calibration_model$calibration_version %in% hybrid_versions) {
        required <- c(required, "he_h2")
    }
    missing <- setdiff(required, names(newdata))
    if (length(missing)) stop("Calibration input is missing: ", paste(missing, collapse = ", "))
    if (calibration_model$calibration_version %in% hybrid_versions) {
        if (identical(calibration_model$calibration_version, "candidate_hybrid_v1")) {
            forward <- predict_candidate_forward_score(
                calibration_model$candidate_model, newdata
            )
        } else {
            forward <- predict_forward_regression(calibration_model$forward_model, newdata)
        }
        forward_h2 <- predict_affine_level_debiasing(
            calibration_model$affine_debiasing, forward$estimate
        )
        estimate_unbounded <- calibration_model$forward_weight * forward_h2 +
            calibration_model$he_weight * newdata$he_h2
        estimate <- pmin(calibration_model$calibration_upper_bound, estimate_unbounded)
        lower <- estimate + calibration_model$residual_quantiles[[1L]]
        upper <- pmin(
            calibration_model$calibration_upper_bound,
            estimate + calibration_model$residual_quantiles[[2L]]
        )
        output <- vector("list", nrow(newdata))
        for (i in seq_len(nrow(newdata))) {
            stratum <- calibration_model$strata[[forward$calibration_stratum[[i]]]]
            raw <- newdata[[raw_metric]][[i]]
            available <- is.finite(raw) && is.finite(estimate[[i]])
            raw_extrapolation <- is.finite(raw) &&
                (raw < stratum$raw_range[[1L]] || raw > stratum$raw_range[[2L]])
            design_extrapolation <-
                forward$calibration_distance[[i]] > calibration_model$max_design_distance
            output[[i]] <- data.frame(
                calibration_stratum = stratum$id,
                calibration_distance = forward$calibration_distance[[i]],
                h2_en_forward = forward_h2[[i]],
                h2_he_weight = calibration_model$he_weight,
                h2_en_calibrated_unbounded = if (available) {
                    estimate_unbounded[[i]]
                } else {
                    NA_real_
                },
                h2_en_calibrated = if (available) estimate[[i]] else NA_real_,
                h2_calibration_lower = if (available) lower[[i]] else NA_real_,
                h2_calibration_upper = if (available) upper[[i]] else NA_real_,
                null_raw_threshold_95 = stratum$null_raw_threshold_95,
                calibration_upper_bound = calibration_model$calibration_upper_bound,
                h2_upper_boundary_hit = if (available) {
                    estimate_unbounded[[i]] > calibration_model$calibration_upper_bound
                } else {
                    NA
                },
                positive_signal = if (is.finite(raw)) {
                    raw > stratum$null_raw_threshold_95
                } else {
                    NA
                },
                calibration_status = if (!is.finite(raw)) {
                    "raw_metric_unavailable"
                } else if (!is.finite(newdata$he_h2[[i]])) {
                    "he_estimate_unavailable"
                } else if (design_extrapolation) {
                    "outside_design_domain"
                } else if (raw_extrapolation) {
                    "raw_metric_extrapolation"
                } else {
                    "within_domain"
                },
                stringsAsFactors = FALSE
            )
        }
        return(do.call(rbind, output))
    }
    output <- vector("list", nrow(newdata))
    for (i in seq_len(nrow(newdata))) {
        distances <- vapply(calibration_model$strata, function(x) {
            calibration_distance(newdata$n[[i]], newdata$num_snps[[i]],
                                 newdata$ld_metric[[i]], x)
        }, numeric(1L))
        selected <- which.min(distances)
        stratum <- calibration_model$strata[[selected]]
        raw <- newdata[[raw_metric]][[i]]
        raw_available <- is.finite(raw)
        estimate <- if (raw_available) {
            stats::approx(
                stratum$curve$raw, stratum$curve$calibrated,
                xout = raw, rule = 2, ties = mean
            )$y
        } else {
            NA_real_
        }
        interval <- estimate + stratum$residual_quantiles
        raw_extrapolation <- raw_available &&
            (raw < stratum$raw_range[[1L]] || raw > stratum$raw_range[[2L]])
        design_extrapolation <- distances[[selected]] > calibration_model$max_design_distance
        output[[i]] <- data.frame(
            calibration_stratum = stratum$id,
            calibration_distance = distances[[selected]],
            h2_en_calibrated = pmin(1, pmax(0, estimate)),
            h2_calibration_lower = pmin(1, pmax(0, interval[[1L]])),
            h2_calibration_upper = pmin(1, pmax(0, interval[[2L]])),
            null_raw_threshold_95 = stratum$null_raw_threshold_95,
            positive_signal = if (raw_available) raw > stratum$null_raw_threshold_95 else NA,
            calibration_status = if (!raw_available) {
                "raw_metric_unavailable"
            } else if (design_extrapolation) {
                "outside_design_domain"
            } else if (raw_extrapolation) {
                "raw_metric_extrapolation"
            } else {
                "within_domain"
            },
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, output)
}

capture_session_info <- function(path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    connection <- file(path, open = "wt")
    on.exit(close(connection), add = TRUE)
    writeLines(c(
        paste("Timestamp:", format(Sys.time(), tz = "UTC")),
        paste("R version:", R.version.string),
        "",
        capture.output(utils::sessionInfo())
    ), connection)
    invisible(path)
}
