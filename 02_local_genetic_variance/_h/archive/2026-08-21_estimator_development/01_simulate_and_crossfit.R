#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    manifest = file.path(dirname(script_path), "..", "_m", "config", "scenarios.tsv"),
    task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = ""),
    output_root = file.path(dirname(script_path), "..", "_m", "raw")
))
task_id <- as_int(cli$task_id, "task_id")
manifest <- read_tsv(cli$manifest)
if (task_id < 1L || task_id > nrow(manifest)) stop("task_id is outside the manifest")
scenario <- manifest[task_id, , drop = FALSE]
if (scenario$scenario_id != task_id) {
    warning("Manifest row and scenario_id differ; using scenario_id from manifest")
}

RNGkind("L'Ecuyer-CMRG")
set.seed(scenario$seed)
simulated <- simulate_locus(
    n = scenario$n,
    p = scenario$num_snps,
    ld_rho = scenario$ld_rho,
    h2 = scenario$true_h2,
    architecture = scenario$architecture
)

fit_error <- NA_character_
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
        fit_error <<- conditionMessage(e)
        NULL
    }
)

if (is.null(fit)) {
    metrics <- data.frame(
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
    metrics <- fit$metrics
}
he_metrics <- haseman_elston(
    simulated$genotype, simulated$phenotype, simulated$covariates
)

result <- cbind(
    scenario[, c(
        "scenario_id", "split", "seed", "n", "num_snps", "ld_rho",
        "architecture", "true_h2", "replicate", "outer_folds",
        "outer_repeats", "inner_folds", "alpha_grid", "lambda_rule",
        "max_features", "raw_metric", "null_alpha", "max_design_distance"
    )],
    metrics[, setdiff(names(metrics), c(
        "n", "num_snps", "outer_folds", "outer_repeats", "inner_folds",
        "max_features"
    )), drop = FALSE]
)
result$realized_h2 <- simulated$realized_h2
result$num_causal <- length(simulated$causal_index)
result$ld_metric <- simulated$ld_metric
result$mean_maf <- mean(colMeans(simulated$genotype) / 2)
result$fit_error <- fit_error
result <- cbind(result, he_metrics)

output <- file.path(
    cli$output_root,
    scenario$split,
    sprintf("scenario-%07d.tsv", scenario$scenario_id)
)
write_tsv(result, output)
cat("Wrote", output, "\n")
