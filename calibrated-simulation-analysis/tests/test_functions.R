#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "..", "_h", "00_functions.R"))

RNGkind("L'Ecuyer-CMRG")
set.seed(1234)
simulation <- simulate_locus(
    n = 54, p = 30, ld_rho = 0.4, h2 = 0.3, architecture = "sparse"
)
stopifnot(abs(simulation$realized_h2 - 0.3) < 1e-8)
stopifnot(is.finite(simulation$ld_metric))

fit_once <- function() crossfit_elastic_net(
    simulation$genotype,
    simulation$phenotype,
    simulation$covariates,
    outer_folds = 3,
    outer_repeats = 1,
    inner_folds = 3,
    alpha_grid = c(0.5, 1),
    max_features = 30,
    seed = 99,
    keep_predictions = TRUE
)
first <- fit_once()
second <- fit_once()
stopifnot(identical(first$metrics, second$metrics))
stopifnot(identical(first$predictions, second$predictions))
stopifnot(all(first$predictions$prediction_repeats == 1L))
stopifnot(nrow(first$folds) == 3L)
stopifnot(is.finite(first$metrics$rho2_oof))
stopifnot(first$metrics$r2_oof <= 1)
he <- haseman_elston(
    simulation$genotype, simulation$phenotype, simulation$covariates
)
stopifnot(he$he_converged)
stopifnot(is.finite(he$he_h2), is.finite(he$he_se), is.finite(he$he_pvalue))
null_cutoff <- finite_sample_upper_threshold(seq_len(30), alpha = 0.05)
stopifnot(null_cutoff$threshold == 30)
stopifnot(null_cutoff$order_index == 30L)
stopifnot(abs(null_cutoff$attainable_alpha - 1 / 31) < 1e-12)
cat("All function tests passed\n")
