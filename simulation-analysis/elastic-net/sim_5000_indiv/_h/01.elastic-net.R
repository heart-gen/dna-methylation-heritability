## --- Load libraries --- ##
library(here)
library(glmnet)
library(bigsnpr)
library(bigstatsr)

## --- Helper Functions --- ##
get_phenotypes <- function(num_samples) {
    base_dir <- here("inputs/simulated-data/_m",
                       paste0("sim_", num_samples, "_indiv"))
    pheno_file <- here(base_dir, "simulated.phen")
    if (!file.exists(pheno_file)) {
        stop("VMR list file not found: ", pheno_file)
    }
    return(read.table(pheno_file, header=TRUE, stringsAsFactors=FALSE))
}

load_genotypes <- function(num_samples) {
    base_dir <- here("inputs/simulated-data/_m",
                     paste0("sim_", num_samples, "_indiv"))
    geno_bed_path <- here(base_dir, "plink_sim/simulated.bed")
    cat("Processing PLINK file:", basename(geno_bed_path), "\n")
                                        # Use tempfile for backingfile to
                                        # avoid conflicts in array jobs
    ##gcc_dir     <- "/projects/b1042/HEART-GeN-Lab/"
    backing_rds <- tempfile(fileext = ".rds")
    rds_path    <- snp_readBed(geno_bed_path,
                               backingfile=sub("\\.rds$", "", backing_rds))
    return(snp_attach(rds_path))
}

get_pheno_loc <- function(num_samples, task_id){
    base_dir <- here("inputs/simulated-data/_m",
                     paste0("sim_", num_samples, "_indiv"))
    mapping_file <- here(base_dir, "snp_phenotype_mapping.tsv")
    mapped_df <- read.table(mapping_file, header=TRUE, stringsAsFactors = FALSE)
    return(mapped_df[(mapped_df["phenotype_id"] == paste0("pheno_", task_id)), ])
}

perform_snp_clumping <- function(G_imputed, info, pheno_scaled) {
                                        # Compute association statistics
    corrs <- big_univLinReg(G_imputed, pheno_scaled)
    stat  <- abs(corrs$estim)
                                        # Perform clumping on sorted data
    ind_keep <- snp_clumping(
        G = G_imputed,
        infos.chr = info$chromosome,
        infos.pos = info$physical.pos,
        S = stat
    )
    return(ind_keep)
}

## --- MAIN SCRIPT --- ##
                                        # Retrieve variables
NUM_SAMPLES <- Sys.getenv("NUM_SAMPLES")
task_id     <- as.integer(Sys.getenv("task_id"))

if (is.na(task_id)) {
    stop("SLURM_ARRAY_TASK_ID is not set or is not a valid integer.")
}
if (NUM_SAMPLES == "") {
    stop("NUM_SAMPLES environment variable is not set.")
}

                                        # Set reproducible seed per task
RNGkind("L'Ecuyer-CMRG")
set.seed(20250525 + task_id)

                                        # Load PHENO data
pheno_list <- get_phenotypes(NUM_SAMPLES)
if (task_id < 1 || task_id+2 > ncol(pheno_list)) {
    stop("SLURM_ARRAY_TASK_ID is out of bounds for the PHENO list.")
}
pheno_entry  <- pheno_list[, task_id+2]

                                        # Load genotypes
bigSNP      <- load_genotypes(NUM_SAMPLES)
bk_file     <- bigSNP$genotypes$backingfile
rds_file    <- bigSNP$genotypes$rds

                                        # Filter data
pheno_locs <- get_pheno_loc(NUM_SAMPLES, task_id)
ind_loc    <- which(bigSNP$map$chromosome == pheno_locs$chrom &
                    bigSNP$map$physical.pos > pheno_locs$start &
                    bigSNP$map$physical.pos < pheno_locs$end)
G_loc       <- bigSNP$genotypes[, ind_loc, drop = FALSE]
map_loc     <- bigSNP$map[ind_loc, ]

                                        # Reconstruct the bigSNP object
bigSNP_loc  <- list(genotypes = G_loc, map = map_loc, fam = bigSNP$fam)
class(bigSNP_loc) <- class(bigSNP)

                                        # Subset data
G         <- bigSNP_loc$genotypes
infos     <- bigSNP_loc$map

                                        # Ensure SNPs are sorted
sorted_idx  <- order(infos$chromosome, infos$physical.pos)
info_sorted <- infos[sorted_idx, ]
G_sorted    <- G[, sorted_idx, drop = FALSE]

                                        # Load Phenotype
pheno_scaled <- scale(pheno_entry)

                                        # Filter out zero-variance SNPs
snp_variances <- big_apply(
    G,
    function(X, ind) {
        apply(X[, ind, drop = FALSE], 2, function(x) var(x, na.rm=TRUE))
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

                                        # Run SNP clumping
clumped_idx <- perform_snp_clumping(G_imputed, infos_filt, pheno_scaled)
G_clumped   <- as_FBM(G_imputed[, clumped_idx])
infos_filt  <- infos_filt[clumped_idx, ]
cat("Number of SNPs before clumping: ", ncol(G_imputed), "\n")
cat("Number of SNPs after clumping: ", ncol(G_clumped), "\n")

if (ncol(G_clumped) == 0) {
    cat("No SNPs left after clumping. Exiting.\n")
    # Write empty results or handle as needed
    quit(save = "no", status = 0)
}

## --- Boosting framework --- ##
cat("Starting boosting framework...\n")
n_iter     <- 100  # Total boosting iterations
batch_size <- min(1000, ncol(G_clumped)) # SNPs per batch

                                        # Initialize
residuals    <- pheno_scaled
h2_estimates <- numeric(n_iter)
accumulated_betas <- FBM(1, ncol(G_clumped), type = "double", init = 0,
                         backingfile = tempfile())

                                        # Boosting loop
for (iter in 1:n_iter) {
    cat("Boosting iteration: ", iter, "\n")
    if (length(residuals) != nrow(G_clumped)){
        stop("Residuals length does not match genotype matrix rows.")
    }
                                        # Top correlated SNPs
    batch_corrs <- big_univLinReg(G_clumped, residuals)
    selected_snps <- order(abs(batch_corrs$estim), decreasing=TRUE)[1:batch_size]

                                        # Fit batch via elastic net
    X_batch <- as_FBM(G_clumped[, selected_snps])
    cv_fit  <- big_spLinReg(X_batch, residuals, alphas = seq(0.05, 1, 0.05),
                            K = 20)
    fit_summary <- summary(cv_fit)

                                        # Check for convergence issues
    if (any(grepl("not converged", fit_summary))) {
        warning("Model did not converge initeration ", iter)
        h2_estimates[iter] <- 0
        next
    }

                                        # Extract kept indicies
    kept_ind <- attr(cv_fit, "ind.col")

    if (length(kept_ind) > 0) {
        batch_pred <- predict(cv_fit, X_batch)
        residuals  <- residuals - batch_pred
        best_betas <- summary(cv_fit, best.only = TRUE)$beta[[1]]
        global_idx <- selected_snps[kept_ind]

        for(i in seq_along(global_idx)){
            idx <- global_idx[i]
            accumulated_betas[1,idx] <- accumulated_betas[1,idx] + best_betas[i]
        }
                                        # Calculate incremental h2
        h2_estimates[iter] <- var(batch_pred)
    } else {
        h2_estimates[iter] <- 0
    }
    cat("Incremental h2 this iteration: ",
        sprintf("%.5f", h2_estimates[iter]), "\n")

                                        # Early stopping condition
    if (iter > 10 && sd(tail(h2_estimates[1:iter], 5)) < 0.0001) {
        cat("Early stopping criterion met at iteration: ", iter, "\n")
        h2_estimates <- h2_estimates[1:iter] # Trim unused part
        break
    }
}

## --- Final Model Refit and Evaluation --- ##
r_squared_cv <- NA
if (ncol(G_clumped) > 0) {
    cat("Refitting final model using Ridge regression...\n")
    final_model <- cv.glmnet(G_clumped[], pheno_scaled, alpha = 0,
                             nfolds = 5, standardize = TRUE)
    lambda_min  <- final_model$lambda.min
    lambda_idx  <- which(abs(final_model$lambda - lambda_min) < 1e-9)
    if (length(lambda_idx) == 0) {
        lambda_idx <- which.min(abs(final_model$lambda - lambda_min))
    }
    pred_cv <- predict(final_model, G_clumped[], s = "lambda.min")
    valid_idx     <- !is.na(pheno_scaled) & !is.na(pred_cv)
    if (sum(valid_idx) > 1) {
        r_squared_cv <- cor(pheno_scaled[valid_idx], pred_cv[valid_idx])^2
    }
    cat(sprintf("Cross-validated R^2 from Ridge regression: %.4f\n",
                r_squared_cv))
} else {
    cat("Skipping Ridge regression as no SNPs are available.\n")
}

## --- Heritability Estimates from Boosting --- ##
final_accumulated_betas <- accumulated_betas[1, ] # Extract vector of betas

clumped_snp_vars <- big_apply(
    G_clumped,
    function(X, ind) {
        apply(X[, ind, drop = FALSE], 2, var)
    },
    a.combine = "c"
)

h2_unscaled <- sum(final_accumulated_betas^2 * clumped_snp_vars)

## --- Save Results --- ##
task_summary_df <- data.frame(
    pheno_id = paste0("pheno_", task_id),
    num_snps = ifelse(exists("G_clumped"), ncol(G_clumped), 0),
    boosting_iterations_performed = ifelse(exists("iter"), iter, 0),
    h2_unscaled = ifelse(exists("h2_unscaled"), h2_unscaled, NA),
    r_squared_cv = ifelse(exists("r_squared_cv"), r_squared_cv, NA)
)

output_df <- data.frame(
    pheno_id = paste0("pheno_", task_id),
    iteration = seq_along(h2_estimates),
    h2_incremental = h2_estimates
)

betas_df <- data.frame(
    pheno_id = paste0("pheno_", task_id),
    snp_id  = infos_filt$marker.ID,
    beta    = final_accumulated_betas
)

dir.create("summary", recursive = TRUE, showWarnings = FALSE)
write.table(task_summary_df,
            file = sprintf(file.path("summary","task_summary_stats_%d.tsv"),
                           task_id), sep = "\t", quote = F, row.names = F)

dir.create("h2", recursive = TRUE, showWarnings = FALSE)
write.table(output_df,
            file = sprintf(file.path("h2", "h2_estimates_%d.tsv"), task_id),
            sep = "\t", quote = FALSE, row.names = FALSE)

dir.create("betas", recursive = TRUE, showWarnings = FALSE)
write.table(betas_df,
            file = sprintf(file.path("betas", "betas_%d.tsv"), task_id),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Total SNP-based h2 (unscaled): %.4f\n", h2_unscaled))
cat(sprintf("Final h2: %.4f\n", sum(h2_estimates)))

                                        # Clean temporary files
if (file.exists(rds_file)) {
    file.remove(bk_file, rds_file)
    cat("Successfully removed the temporary files.\n")
}

## Reproducibility
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
