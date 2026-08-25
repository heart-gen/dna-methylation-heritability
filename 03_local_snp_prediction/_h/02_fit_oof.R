#!/usr/bin/env Rscript
#### 03_local_snp_prediction -- end-to-end out-of-fold fitting for one chunk ####
##
## Usage (inside an array task):
##   Rscript _h/02_fit_oof.R --run-id lsp-AA-caudate-20260817 --task-ids 1,2,3
##
## This is the script that defect E1 exists to prevent. The legacy code reported
##
##     cor(pheno, predict(cv.glmnet(G_clumped, pheno)))^2
##
## where BOTH the clumping and the lambda choice had seen every donor. Reported
## ~0.85; honest held-out ~0.01.
##
## The rule here is mechanical and absolute: inside an outer fold, the held-out
## donors contribute NOTHING. Not to covariate residualization, not to phenotype
## centering or scaling, not to MAF or missingness or zero-variance filtering,
## not to mean imputation, not to genotype scaling, not to the cis screen, not
## to lambda or alpha selection. Every quantity derived from data is computed on
## the training donors and then APPLIED to the held-out ones. Any new
## preprocessing step added to this file must follow the same
## fit-on-train/apply-to-test shape.
##
## A VMR whose fold-internal screen fails does NOT get dropped: its held-out
## donors receive the prespecified null prediction (the training mean). Dropping
## failed folds would condition the metric on the outcome and reintroduce the
## optimism the module exists to remove.

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(glmnet)
})

MODULE <- "03_local_snp_prediction"
H_DIR <- file.path(repo_root(), MODULE, "_h")

## adjust_phenotype_in_fold() is Module 02's fold-internal residualizer. 03 uses
## the SAME function so that the phenotype 03 predicts is the phenotype 02
## estimated variance on, adjusted the same way. haseman_elston() comes along
## for the screen's equivalence test.
source(file.path(repo_root(), "02_local_genetic_variance", "_h", "00_functions.R"))
source(file.path(H_DIR, "he_permutation_screen.R"))

opts <- parse_v2_args(require = c("run_id", "task_ids"))
allow_unlocked <- isTRUE(opts$allow_unlocked)

run_dir <- file.path(repo_root(), MODULE, "_m", "runs", opts$run_id)
if (!dir.exists(run_dir)) stop("No such run: ", run_dir)

manifest <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
mval <- function(f) {
    v <- manifest$value[manifest$field == f]
    if (length(v) == 0) NA_character_ else v[1]
}
cohort <- mval("cohort"); region <- mval("region")

prediction <- load_config("prediction")
thresholds <- load_config("thresholds")
assert_locked(list(prediction = prediction, thresholds = thresholds),
              allow_unlocked = allow_unlocked)

locked_value <- function(cfg, key, smoke_default) {
    v <- tryCatch(config_get(cfg, key), error = function(e) NULL)
    if (!is.null(v)) return(v)
    if (!allow_unlocked) {
        stop("config/prediction.yml key '", key, "' is null and this is a ",
             "production run (AGENTS.md 12/14).", call. = FALSE)
    }
    smoke_default
}

N_INNER     <- as.integer(locked_value(prediction, "folds.inner", 5L))
ALPHA_GRID  <- as.numeric(locked_value(prediction, "elastic_net.alpha_grid",
                                       c(0.1, 0.5, 0.9, 1.0)))
LAMBDA_RULE <- as.character(locked_value(prediction, "elastic_net.lambda_rule",
                                         "lambda.min"))
SCREEN_METHOD <- as.character(locked_value(prediction, "screen.method",
                                           "he_permutation"))
SCREEN_ALPHA  <- as.numeric(locked_value(prediction, "screen.alpha", 0.05))
N_PERM        <- as.integer(locked_value(prediction, "screen.n_permutations", 1000L))
## The cis QC constants are the SAME keys Module 02 used (config/thresholds.yml
## `cis`), not a 03-specific block: if 03 screened a different variant set than
## 02, the prediction endpoint and the variance endpoint would not be describing
## the same locus. The difference is only WHEN they are applied -- 02 filters on
## all donors, 03 refilters inside each outer training fold.
MAF_MIN       <- as.numeric(locked_value(thresholds, "cis.maf_min", 0.05))
MISS_MAX      <- as.numeric(locked_value(thresholds, "cis.missingness_max", 0.05))
MIN_CIS_VARIANTS <- as.integer(locked_value(thresholds, "cis.min_cis_variants", 100L))

## A smoke run may shrink the permutation count; a production run may not.
if (allow_unlocked) {
    N_PERM <- as.integer(Sys.getenv("LSP_SMOKE_N_PERM", N_PERM))
}

task_ids <- as.integer(strsplit(opts$task_ids, ",", fixed = TRUE)[[1]])
if (any(is.na(task_ids)) || any(task_ids < 1)) {
    stop("--task-ids must be a comma-separated list of positive integers")
}

tasks <- fread(file.path(run_dir, "task-manifest.tsv"))
folds <- fread(file.path(run_dir, "donor-folds.tsv"), colClasses = list(character = "donor"))

out_dir <- file.path(run_dir, "results", "oof")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fail_dir <- file.path(run_dir, "results", "failures")
dir.create(fail_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------------ inputs
## Module 01 runs have no `results/` directory: genotypes are per-chromosome
## PLINK1 sets, phenotypes are one .phen per VMR, covariates are per-chromosome
## .covar/.qcovar. 00_shared/locus_io.R is the single reader for that layout and
## is shared verbatim with Module 02 Stage 01, so 02 and 03 provably score the
## same variants, the same donors, and the same covariate design.
vmr_run <- mval("upstream_vmr_catalog_run_id")
if (is.na(vmr_run) || !nzchar(vmr_run)) {
    stop("Run manifest carries no upstream_vmr_catalog_run_id; rerun 00_new_run.R")
}
vmr_run_dir <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs", vmr_run)
if (!dir.exists(vmr_run_dir)) stop("Upstream 01 run not found: ", vmr_run_dir)

#' Fit-on-train / apply-to-test genotype preprocessing.
#'
#' Returns NULL when nothing survives the TRAINING-set filters, which is a
#' screen failure, not an error.
prep_genotypes <- function(g_train, g_test) {
    ## Missingness and MAF are computed on training donors only.
    miss <- colMeans(is.na(g_train))
    keep <- miss <= MISS_MAX
    if (!any(keep)) return(NULL)
    g_train <- g_train[, keep, drop = FALSE]
    g_test  <- g_test[, keep, drop = FALSE]

    train_mean <- colMeans(g_train, na.rm = TRUE)
    maf <- pmin(train_mean, 2 - train_mean) / 2
    keep <- !is.na(maf) & maf >= MAF_MIN
    if (!any(keep)) return(NULL)
    g_train <- g_train[, keep, drop = FALSE]
    g_test  <- g_test[, keep, drop = FALSE]
    train_mean <- train_mean[keep]

    ## Mean-impute BOTH matrices with the training means. Imputing the test
    ## matrix from its own column means is a classic silent leak.
    for (j in seq_len(ncol(g_train))) {
        na_tr <- is.na(g_train[, j]); if (any(na_tr)) g_train[na_tr, j] <- train_mean[j]
        na_te <- is.na(g_test[, j]);  if (any(na_te)) g_test[na_te, j]  <- train_mean[j]
    }

    train_sd <- apply(g_train, 2, stats::sd)
    keep <- train_sd > 0
    if (!any(keep)) return(NULL)
    g_train <- g_train[, keep, drop = FALSE]
    g_test  <- g_test[, keep, drop = FALSE]
    train_mean <- train_mean[keep]; train_sd <- train_sd[keep]

    list(
        train = scale(g_train, center = train_mean, scale = train_sd),
        test  = scale(g_test,  center = train_mean, scale = train_sd),
        n_variants = ncol(g_train)
    )
}

#' Fold-internal cis screen. Computed on TRAINING donors only.
#'
#' `he_permutation` is the locked production method (config/prediction.yml).
#' The bare `haseman_elston()` p-value is NOT offered as an option: AGENTS.md
#' 7.3 forbids it uncalibrated, and an unreachable branch is safer than a
#' reachable wrong one.
screen_locus <- function(y_train, x_train, seed) {
    if (SCREEN_METHOD %in% c("none", "")) {
        return(list(pass = TRUE, p = NA_real_, statistic = NA_real_))
    }
    if (SCREEN_METHOD == "he_permutation") {
        prep <- he_screen_prepare(x_train)
        return(he_permutation_screen(prep, y_train, n_perm = N_PERM,
                                     alpha = SCREEN_ALPHA, seed = seed))
    }
    if (SCREEN_METHOD == "marginal_min_p") {
        p <- apply(x_train, 2, function(v) {
            fit <- stats::lm(y_train ~ v)
            s <- summary(fit)$coefficients
            if (nrow(s) < 2) NA_real_ else s[2, 4]
        })
        pmin_adj <- min(stats::p.adjust(p, method = "bonferroni"), na.rm = TRUE)
        return(list(pass = is.finite(pmin_adj) && pmin_adj <= SCREEN_ALPHA,
                    p = pmin_adj, statistic = NA_real_))
    }
    stop("Unsupported screen.method '", SCREEN_METHOD, "'. Implement it here ",
         "explicitly rather than falling through to an unscreened fit.")
}

## -------------------------------------------------------------------- loop
for (tid in task_ids) {
    row <- tasks[tasks$task_id == tid]
    if (nrow(row) != 1) stop("task_id ", tid, " not in the task manifest")
    vmr_id <- row$vmr_id[1]
    out_f <- file.path(out_dir, paste0(vmr_id, ".tsv"))
    if (file.exists(out_f)) next          # idempotent restart

    fit_one <- function() {
        locus <- load_observed_locus(
            task = list(chrom = row$chrom[1], start = row$start[1],
                        end = row$end[1]),
            cohort = cohort, vmr_run_dir = vmr_run_dir,
            min_cis_variants = MIN_CIS_VARIANTS, backing_tag = "lsp",
            ## Every data-dependent genotype filter happens inside the fold.
            apply_snp_qc = FALSE
        )
        if (!identical(locus$status, "ok")) {
            ## Excluded or QC-failed upstream is a documented outcome, not a
            ## computational failure. Record it so 03's reconciliation can
            ## account for the task without counting it as a crash.
            writeLines(paste0("status=", locus$status, " reason=", locus$reason),
                       file.path(out_dir, paste0(vmr_id, ".skip")))
            return(paste0(locus$status, " (", locus$reason, ")"))
        }

        donor_key <- paste(locus$metadata$FID, locus$metadata$IID, sep = "::")
        y_all <- stats::setNames(as.numeric(locus$y), donor_key)
        g_all <- locus$genotype
        rownames(g_all) <- donor_key
        covar_all <- as.matrix(locus$covariates)
        rownames(covar_all) <- donor_key

        shared <- intersect(donor_key, unique(folds$donor))
        if (length(shared) < 20) stop("only ", length(shared), " usable donors")
        y_all <- y_all[shared]
        g_all <- g_all[shared, , drop = FALSE]
        covar_all <- covar_all[shared, , drop = FALSE]

        preds <- rbindlist(lapply(split(folds, folds$repeat_i), function(fr) {
            rbindlist(lapply(sort(unique(fr$outer_fold)), function(k) {
                test_ids  <- intersect(fr$donor[fr$outer_fold == k], shared)
                train_ids <- setdiff(shared, test_ids)
                if (length(test_ids) == 0 || length(train_ids) < 10) return(NULL)

                ## Covariate residualization AND centering/scaling, both learned
                ## on training donors only and applied to the held-out ones.
                ## Without this, age, sex and diagnosis are uncontrolled and 03's
                ## R2 is not comparable to 02's variance estimate.
                adj <- tryCatch(
                    adjust_phenotype_in_fold(
                        y_train = y_all[train_ids], y_test = y_all[test_ids],
                        covar_train = covar_all[train_ids, , drop = FALSE],
                        covar_test = covar_all[test_ids, , drop = FALSE]),
                    error = function(e) NULL)
                if (is.null(adj)) return(NULL)
                y_tr <- adj$train; y_te <- adj$test

                gp <- prep_genotypes(g_all[train_ids, , drop = FALSE],
                                     g_all[test_ids, , drop = FALSE])

                ## The prespecified null prediction: the training mean, which is
                ## 0 on the training-residualized, training-standardized scale.
                null_row <- data.table(
                    repeat_i = fr$repeat_i[1], outer_fold = k,
                    donor = test_ids, y_obs = y_te, y_pred = 0,
                    screened_in = FALSE, n_variants = 0L,
                    screen_p = NA_real_, screen_stat = NA_real_,
                    alpha = NA_real_, lambda = NA_real_)

                if (is.null(gp)) return(null_row)

                fold_seed <- seed_for(opts$run_id, region = region, task = vmr_id,
                                      repeat_i = fr$repeat_i[1], fold = k)
                sc <- screen_locus(y_tr, gp$train, seed = fold_seed)
                if (!isTRUE(sc$pass)) {
                    null_row[, screen_p := sc$p]
                    null_row[, screen_stat := sc$statistic]
                    null_row[, n_variants := gp$n_variants]
                    return(null_row)
                }

                ## Inner CV picks alpha and lambda on training donors only. The
                ## same inner fold ids are reused across the alpha grid so the
                ## comparison between alphas is paired.
                set.seed(fold_seed)
                inner_id <- sample(rep(seq_len(N_INNER), length.out = length(train_ids)))

                fits <- lapply(ALPHA_GRID, function(a) {
                    tryCatch(cv.glmnet(gp$train, y_tr, alpha = a,
                                       foldid = inner_id, standardize = FALSE),
                             error = function(e) NULL)
                })
                ok <- !vapply(fits, is.null, logical(1))
                if (!any(ok)) return(null_row)
                fits <- fits[ok]; alphas <- ALPHA_GRID[ok]

                cvm <- vapply(seq_along(fits), function(i) {
                    f <- fits[[i]]
                    f$cvm[match(f[[LAMBDA_RULE]], f$lambda)]
                }, numeric(1))
                best <- which.min(cvm)
                f <- fits[[best]]
                lam <- f[[LAMBDA_RULE]]

                yhat <- as.numeric(predict(f, newx = gp$test, s = lam))
                data.table(repeat_i = fr$repeat_i[1], outer_fold = k,
                           donor = test_ids, y_obs = y_te, y_pred = yhat,
                           screened_in = TRUE, n_variants = gp$n_variants,
                           screen_p = sc$p, screen_stat = sc$statistic,
                           alpha = alphas[best], lambda = lam)
            }))
        }))

        if (is.null(preds) || nrow(preds) == 0) stop("no fold produced predictions")
        preds[, vmr_id := vmr_id]
        write_atomic(preds, out_f)
        "ok"
    }

    result <- tryCatch(fit_one(), error = function(e) {
        writeLines(conditionMessage(e), file.path(fail_dir, paste0(vmr_id, ".txt")))
        "failed"
    })

    message("[03] task ", tid, " (", vmr_id, "): ", result)
}
