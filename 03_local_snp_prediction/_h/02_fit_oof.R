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
## donors contribute NOTHING. Not to phenotype centering or scaling, not to MAF
## or missingness or zero-variance filtering, not to mean imputation, not to
## genotype scaling, not to the cis screen, not to lambda or alpha selection.
## Every quantity derived from data is computed on the training donors and then
## APPLIED to the held-out ones. Any new preprocessing step added to this file
## must follow the same fit-on-train/apply-to-test shape.
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
                                       c(0.1, 0.5, 0.9)))
LAMBDA_RULE <- as.character(locked_value(prediction, "elastic_net.lambda_rule",
                                         "lambda.min"))
SCREEN_METHOD <- as.character(locked_value(prediction, "screen.method", "none"))
SCREEN_ALPHA  <- as.numeric(locked_value(prediction, "screen.multiple_testing", 1))
MAF_MIN       <- as.numeric(locked_value(thresholds, "prediction.maf_min", 0.01))
MISS_MAX      <- as.numeric(locked_value(thresholds, "prediction.missingness_max", 0.05))

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
vmr_run <- mval("upstream_vmr_catalog_run_id")
pheno <- fread(file.path(repo_root(), "01_vmr_catalog", "_m", "runs", vmr_run,
                         "results", "vmr_meth.phen"))
pheno_donors <- as.character(pheno[[2]])

#' Read the cis genotype dosage matrix for one VMR.
#'
#' 01_vmr_catalog's genotype-extraction step already wrote a per-VMR plink set
#' for the locked cis window; 03 reuses it rather than re-extracting, so 02 and
#' 03 are guaranteed to be scoring the same variants. plink2 comes from the opt
#' tree, never the module system (defect V12).
read_cis_dosages <- function(vmr_id) {
    prefix <- file.path(repo_root(), "01_vmr_catalog", "_m", "runs", vmr_run,
                        "results", "plink_format", vmr_id)
    if (!file.exists(paste0(prefix, ".pgen"))) {
        return(NULL)
    }
    raw <- tempfile(fileext = ".raw")
    on.exit(unlink(c(raw, sub("\\.raw$", ".log", raw))), add = TRUE)
    plink2 <- Sys.getenv("PLINK2", "/projects/p32505/opt/bin/plink2")
    status <- system2(plink2, c("--pfile", shQuote(prefix),
                                "--export", "A",
                                "--threads", Sys.getenv("V2_THREADS", "1"),
                                "--out", shQuote(sub("\\.raw$", "", raw))),
                      stdout = FALSE, stderr = FALSE)
    if (status != 0 || !file.exists(raw)) return(NULL)
    g <- fread(raw)
    ids <- as.character(g$IID)
    mat <- as.matrix(g[, -(1:6)])
    rownames(mat) <- ids
    mat
}

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
#' NOTE: config/prediction.yml deliberately ships screen.method: null. The
#' legacy ordinary-OLS Haseman-Elston p-value must not become the production
#' screen without donor-robust or permutation calibration (README, AGENTS.md
#' 7.2), so "none" is the only method wired up until the PI locks one.
screen_locus <- function(y_train, x_train) {
    if (SCREEN_METHOD %in% c("none", "")) return(list(pass = TRUE, p = NA_real_))
    if (SCREEN_METHOD == "marginal_min_p") {
        p <- apply(x_train, 2, function(v) {
            fit <- stats::lm(y_train ~ v)
            s <- summary(fit)$coefficients
            if (nrow(s) < 2) NA_real_ else s[2, 4]
        })
        pmin_adj <- min(stats::p.adjust(p, method = "bonferroni"), na.rm = TRUE)
        return(list(pass = is.finite(pmin_adj) && pmin_adj <= SCREEN_ALPHA,
                    p = pmin_adj))
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

    result <- tryCatch({
        col <- match(vmr_id, names(pheno))
        if (is.na(col)) stop("VMR ", vmr_id, " absent from the phenotype matrix")
        y_all <- stats::setNames(as.numeric(pheno[[col]]), pheno_donors)

        g_all <- read_cis_dosages(vmr_id)
        if (is.null(g_all)) stop("no cis genotypes available for ", vmr_id)

        ## align_by_id is the V1 fix: never assume two sources share row order.
        shared <- intersect(names(y_all), rownames(g_all))
        shared <- intersect(shared, unique(folds$donor))
        if (length(shared) < 20) stop("only ", length(shared), " usable donors")
        y_all <- y_all[shared]
        g_all <- g_all[shared, , drop = FALSE]

        preds <- rbindlist(lapply(split(folds, folds$repeat_i), function(fr) {
            rbindlist(lapply(sort(unique(fr$outer_fold)), function(k) {
                test_ids  <- intersect(fr$donor[fr$outer_fold == k], shared)
                train_ids <- setdiff(shared, test_ids)
                if (length(test_ids) == 0 || length(train_ids) < 10) return(NULL)

                ## Phenotype centering/scaling from TRAINING donors only.
                y_tr_raw <- y_all[train_ids]
                mu <- mean(y_tr_raw); sdev <- stats::sd(y_tr_raw)
                if (!is.finite(sdev) || sdev == 0) return(NULL)
                y_tr <- (y_tr_raw - mu) / sdev
                y_te <- (y_all[test_ids] - mu) / sdev

                gp <- prep_genotypes(g_all[train_ids, , drop = FALSE],
                                     g_all[test_ids, , drop = FALSE])

                ## The prespecified null prediction: the training mean, which is
                ## 0 on the training-standardized scale.
                null_row <- data.table(
                    repeat_i = fr$repeat_i[1], outer_fold = k,
                    donor = test_ids, y_obs = y_te, y_pred = 0,
                    screened_in = FALSE, n_variants = 0L,
                    screen_p = NA_real_, alpha = NA_real_, lambda = NA_real_)

                if (is.null(gp)) return(null_row)

                sc <- screen_locus(y_tr, gp$train)
                if (!isTRUE(sc$pass)) {
                    null_row[, screen_p := sc$p]
                    null_row[, n_variants := gp$n_variants]
                    return(null_row)
                }

                ## Inner CV picks alpha and lambda on training donors only. The
                ## same inner fold ids are reused across the alpha grid so the
                ## comparison between alphas is paired.
                set.seed(seed_for(opts$run_id, region = region, task = vmr_id,
                                  repeat_i = fr$repeat_i[1], fold = k))
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
                           screen_p = sc$p, alpha = alphas[best], lambda = lam)
            }))
        }))

        if (is.null(preds) || nrow(preds) == 0) stop("no fold produced predictions")
        preds[, vmr_id := vmr_id]
        write_atomic(preds, out_f)
        "ok"
    }, error = function(e) {
        writeLines(conditionMessage(e), file.path(fail_dir, paste0(vmr_id, ".txt")))
        "failed"
    })

    message("[03] task ", tid, " (", vmr_id, "): ", result)
}
