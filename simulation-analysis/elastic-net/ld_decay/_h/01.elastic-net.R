## --- Load libraries --- ##
library(here)
library(glmnet)
library(bigsnpr)
library(bigstatsr)

## --- Helper Functions --- ##
get_sim_input_dir <- function() {
    run_name <- Sys.getenv("RUN_NAME")
    sim_input_dir <- here("inputs/simulated-data/_m", run_name)
    if (run_name == "") {
        stop("RUN_NAME environment variable is not set.")
    }
    return(sim_input_dir)
}

get_phenotypes <- function(sim_input_dir) {
    pheno_file <- file.path(sim_input_dir, "simulated.phen")
    if (!file.exists(pheno_file)) {
        stop("Phenotype file not found: ", pheno_file)
    }
    return(read.table(pheno_file, header=TRUE, stringsAsFactors=FALSE))
}

load_genotypes <- function(sim_input_dir) {
    geno_bed_path <- file.path(sim_input_dir, "plink_sim", "simulated.bed")
    if (!file.exists(geno_bed_path)) {
        stop("Genotype BED file not found: ", geno_bed_path)
    }
    cat("Processing PLINK file:", basename(geno_bed_path), "\n")
                                        # Use tempfile for backingfile to
                                        # avoid conflicts in array jobs
    backing_rds <- tempfile(fileext = ".rds")
    rds_path    <- snp_readBed(geno_bed_path,
                               backingfile=sub("\\.rds$", "", backing_rds))
    return(snp_attach(rds_path))
}

get_pheno_loc <- function(sim_input_dir, task_id){
    mapping_file <- file.path(sim_input_dir, "snp_phenotype_mapping.tsv")
    if (!file.exists(mapping_file)) {
        stop("Phenotype mapping file not found: ", mapping_file)
    }
    mapped_df <- read.table(mapping_file, header=TRUE, stringsAsFactors = FALSE)
    return(mapped_df[(mapped_df["phenotype_id"] == paste0("pheno_", task_id)), ])
}

## --- Method implementations --- ##

run_boosting_hybrid <- function(G_imputed, pheno_scaled) {
    ## Boosting loop: screens active SNP set from all QC-passed SNPs.
    ## h2 and r2 are estimated from a joint ridge refit on the active set,
    ## NOT from summed boosting increments.
    cat("METHOD: boosting_hybrid\n")
    n_iter     <- 100
    batch_size <- min(1000, ncol(G_imputed))

    residuals         <- pheno_scaled
    h2_iters          <- numeric(n_iter)
    accumulated_betas <- FBM(1, ncol(G_imputed), type = "double", init = 0,
                             backingfile = tempfile())
    iters_done <- 0

    for (iter in 1:n_iter) {
        cat("Boosting iteration:", iter, "\n")
        if (length(residuals) != nrow(G_imputed)) {
            stop("Residuals length does not match genotype matrix rows.")
        }
                                        # Select top correlated SNPs
        batch_corrs   <- big_univLinReg(G_imputed, residuals)
        selected_snps <- order(abs(batch_corrs$estim), decreasing = TRUE)[1:batch_size]

                                        # Fit elastic net on batch
        X_batch    <- as_FBM(G_imputed[, selected_snps])
        cv_fit     <- big_spLinReg(X_batch, residuals,
                                   alphas = seq(0.05, 1, 0.05), K = 5)
        kept_batch <- attr(cv_fit, "ind.col")
        iters_done <- iter

        if (length(kept_batch) > 0) {
            batch_pred <- predict(cv_fit, X_batch)
            residuals  <- residuals - batch_pred
            h2_iters[iter] <- var(batch_pred)
            betas_b    <- summary(cv_fit, best.only = TRUE)$beta[[1]]
            global_idx <- selected_snps[kept_batch]
            for (i in seq_along(global_idx)) {
                accumulated_betas[1, global_idx[i]] <-
                    accumulated_betas[1, global_idx[i]] + betas_b[i]
            }
        } else {
            h2_iters[iter] <- 0
        }
        cat("Incremental h2 this iteration:",
            sprintf("%.5f", h2_iters[iter]), "\n")

                                        # Early stopping
        if (iter > 10 &&
            sd(tail(h2_iters[1:iter], 5)) < 0.0001) {
            cat("Early stopping criterion met at iteration:", iter, "\n")
            h2_iters <- h2_iters[1:iter]
            break
        }
    }

    active_idx <- which(accumulated_betas[1, ] != 0)
    cat("Active SNPs from boosting:", length(active_idx), "\n")

                                        # Joint ridge refit on active set
    r_squared_cv <- NA
    final_betas  <- rep(0, ncol(G_imputed))

    if (length(active_idx) > 0) {
        cat("Refitting joint ridge on active set...\n")
        X_active    <- G_imputed[, active_idx, drop = FALSE]
        ridge_model <- cv.glmnet(X_active, pheno_scaled, alpha = 0,
                                 nfolds = 5, standardize = TRUE)
        final_betas[active_idx] <- as.numeric(
            coef(ridge_model, s = "lambda.min")[-1]
        )
        pred_cv   <- predict(ridge_model, X_active, s = "lambda.min")
        valid_idx <- !is.na(pheno_scaled) & !is.na(pred_cv)
        if (sum(valid_idx) > 1) {
            r_squared_cv <- cor(pheno_scaled[valid_idx], pred_cv[valid_idx])^2
        }
        cat(sprintf("Joint-ridge R^2: %.4f\n", r_squared_cv))
    } else {
        cat("No active SNPs — skipping ridge refit.\n")
    }

    list(
        final_betas          = final_betas,
        active_idx           = active_idx,
        r_squared_cv         = r_squared_cv,
        iters_done           = iters_done,
        h2_iters             = h2_iters
    )
}

run_joint_ridge <- function(G_imputed, pheno_scaled) {
    ## Single big_spLinReg call with ridge-dominant alphas.
    ## Screening is applied within each training fold by bigstatsr.
    cat("METHOD: joint_ridge\n")

    cat("Fitting joint elastic net (alphas 1e-4 to 0.2, fold-internal screening)...\n")
    cv_fit <- big_spLinReg(G_imputed, pheno_scaled,
                           alphas = c(1e-4, 0.05, 0.1, 0.15, 0.2), K = 5)

    kept_ind       <- attr(cv_fit, "ind.col")
    final_betas    <- rep(0, ncol(G_imputed))
    if (length(kept_ind) > 0) {
        final_betas[kept_ind] <- summary(cv_fit, best.only = TRUE)$beta[[1]]
    }
    active_idx <- kept_ind
    cat("SNPs with non-zero betas:", length(active_idx), "\n")

                                        # Ridge refit for r2
    r_squared_cv <- NA
    if (length(active_idx) > 0) {
        cat("Computing cross-validated R^2 via ridge regression...\n")
        X_kept      <- G_imputed[, active_idx, drop = FALSE]
        ridge_model <- cv.glmnet(X_kept, pheno_scaled, alpha = 0,
                                 nfolds = 5, standardize = TRUE)
        pred_cv   <- predict(ridge_model, X_kept, s = "lambda.min")
        valid_idx <- !is.na(pheno_scaled) & !is.na(pred_cv)
        if (sum(valid_idx) > 1) {
            r_squared_cv <- cor(pheno_scaled[valid_idx], pred_cv[valid_idx])^2
        }
        cat(sprintf("Joint-ridge R^2: %.4f\n", r_squared_cv))
    } else {
        cat("No SNPs selected — skipping r^2 computation.\n")
    }

    list(
        final_betas          = final_betas,
        active_idx           = active_idx,
        r_squared_cv         = r_squared_cv,
        iters_done           = 1L,
        h2_iters             = NA_real_
    )
}

## --- MAIN SCRIPT --- ##
                                        # Retrieve variables
METHOD      <- Sys.getenv("METHOD", unset = "boosting_hybrid")
NUM_SAMPLES <- Sys.getenv("NUM_SAMPLES")
task_id     <- as.integer(Sys.getenv("task_id"))
SIM_INPUT_DIR <- get_sim_input_dir()

if (is.na(task_id)) {
    stop("task_id is not set or is not a valid integer.")
}
if (NUM_SAMPLES == "") {
    stop("NUM_SAMPLES environment variable is not set.")
}
if (!METHOD %in% c("boosting_hybrid", "joint_ridge")) {
    stop("METHOD must be 'boosting_hybrid' or 'joint_ridge', got: ", METHOD)
}

cat("Using simulation input directory:", SIM_INPUT_DIR, "\n")
cat("METHOD:", METHOD, "\n")

                                        # Set reproducible seed per task
RNGkind("L'Ecuyer-CMRG")
set.seed(20250525 + task_id)

                                        # Load PHENO data
pheno_list <- get_phenotypes(SIM_INPUT_DIR)
if (task_id < 1 || task_id + 2 > ncol(pheno_list)) {
    stop("task_id is out of bounds for the PHENO list.")
}
pheno_entry  <- pheno_list[, task_id + 2]

                                        # Load genotypes
bigSNP  <- load_genotypes(SIM_INPUT_DIR)
bk_file <- bigSNP$genotypes$backingfile
rds_file <- bigSNP$genotypes$rds

                                        # Filter to phenotype locus
pheno_locs <- get_pheno_loc(SIM_INPUT_DIR, task_id)
ind_loc    <- which(bigSNP$map$chromosome == pheno_locs$chrom &
                    bigSNP$map$physical.pos > pheno_locs$start &
                    bigSNP$map$physical.pos < pheno_locs$end)
G_loc      <- bigSNP$genotypes[, ind_loc, drop = FALSE]
map_loc    <- bigSNP$map[ind_loc, ]

                                        # Reconstruct localised bigSNP object
bigSNP_loc        <- list(genotypes = G_loc, map = map_loc, fam = bigSNP$fam)
class(bigSNP_loc) <- class(bigSNP)

                                        # Subset data
G     <- bigSNP_loc$genotypes
infos <- bigSNP_loc$map

                                        # Sort SNPs
sorted_idx  <- order(infos$chromosome, infos$physical.pos)
info_sorted <- infos[sorted_idx, ]
G_sorted    <- G[, sorted_idx, drop = FALSE]

                                        # Scale phenotype
pheno_scaled <- scale(pheno_entry)

                                        # Filter zero-variance SNPs
snp_variances <- big_apply(
    G,
    function(X, ind) {
        apply(X[, ind, drop = FALSE], 2, function(x) var(x, na.rm = TRUE))
    },
    a.combine = "c"
)
keep_idx   <- which(snp_variances > 1e-6)
infos_filt <- info_sorted[keep_idx, ]
G_temp     <- G_sorted[, keep_idx]
G_filtered <- FBM.code256(
    nrow = nrow(G_temp), ncol = ncol(G_temp),
    code = bigSNP$genotypes$code256,
    backingfile = tempfile()
)
G_filtered[] <- G_temp[]
rm(G_loc, map_loc, bigSNP, G_temp, G, infos, G_sorted, info_sorted)

                                        # Impute missing values
cat("Imputing missing genotypes using mode...\n")
G_imputed <- snp_fastImputeSimple(G_filtered, method = "mode")

cat("Number of SNPs after QC:", ncol(G_imputed), "\n")

if (ncol(G_imputed) == 0) {
    cat("No SNPs available after QC. Exiting.\n")
    quit(save = "no", status = 0)
}

## --- Dispatch to selected method --- ##
result <- if (METHOD == "boosting_hybrid") {
    run_boosting_hybrid(G_imputed, pheno_scaled)
} else {
    run_joint_ridge(G_imputed, pheno_scaled)
}

final_betas  <- result$final_betas
active_idx   <- result$active_idx
r_squared_cv <- result$r_squared_cv
iters_done   <- result$iters_done
h2_iters     <- result$h2_iters

## --- Heritability Estimation --- ##
## h2 is always computed from joint-ridge betas, not boosting increments
snp_vars <- big_apply(
    G_imputed,
    function(X, ind) apply(X[, ind, drop = FALSE], 2, var),
    a.combine = "c"
)
h2_unscaled <- sum(final_betas^2 * snp_vars)

## --- Save Results --- ##
task_summary_df <- data.frame(
    pheno_id = paste0("pheno_", task_id),
    num_snps = length(active_idx),
    boosting_iterations_performed = iters_done,
    h2_unscaled  = h2_unscaled,
    r_squared_cv = r_squared_cv
)

                                        # h2 per iteration (NA row for joint_ridge)
if (all(is.na(h2_iters))) {
    output_df <- data.frame(
        pheno_id       = paste0("pheno_", task_id),
        iteration      = 1L,
        h2_incremental = h2_unscaled
    )
} else {
    output_df <- data.frame(
        pheno_id       = paste0("pheno_", task_id),
        iteration      = seq_along(h2_iters),
        h2_incremental = h2_iters
    )
}

                                        # Betas: only non-zero SNPs
betas_df <- data.frame(
    pheno_id = paste0("pheno_", task_id),
    snp_id   = infos_filt$marker.ID[active_idx],
    beta     = final_betas[active_idx]
)

                                        # Use method-specific subdirectories so
                                        # boosting_hybrid and joint_ridge runs
                                        # can proceed concurrently without
                                        # overwriting each other's task files.
summary_dir <- paste0("summary_", METHOD)
h2_dir      <- paste0("h2_",      METHOD)
betas_dir   <- paste0("betas_",   METHOD)

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
write.table(task_summary_df,
            file = sprintf(file.path(summary_dir, "task_summary_stats_%d.tsv"),
                           task_id),
            sep = "\t", quote = FALSE, row.names = FALSE)

dir.create(h2_dir, recursive = TRUE, showWarnings = FALSE)
write.table(output_df,
            file = sprintf(file.path(h2_dir, "h2_estimates_%d.tsv"), task_id),
            sep = "\t", quote = FALSE, row.names = FALSE)

dir.create(betas_dir, recursive = TRUE, showWarnings = FALSE)
write.table(betas_df,
            file = sprintf(file.path(betas_dir, "betas_%d.tsv"), task_id),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Total SNP-based h2 (unscaled): %.4f\n", h2_unscaled))

                                        # Clean temporary files
if (file.exists(rds_file)) {
    file.remove(bk_file, rds_file)
    cat("Successfully removed temporary files.\n")
}

## --- Reproducibility --- ##
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
