## --- Load libraries --- ##
library(here)
library(glmnet)
library(bigsnpr)
library(bigstatsr)

## --- Helper Functions --- ##
get_error_list <- function(error_file_path = "../_h/snp-error-window.tsv") {
                                        # Read error regions
    if (file.exists(error_file_path)) {
        error_regions <- read.table(error_file_path,header=TRUE,sep="\t") |>
            dplyr::mutate(Chrom=as.numeric(gsub("chr", "", Chr)),
                          Start=as.numeric(Start), End=as.numeric(End)) |>
            dplyr::select(Chrom, Start, End)
    } else {
        message("Warning: Error regions file not found at: ",
                error_file_path)
        error_regions <- data.frame(Chrom = numeric(), Start = numeric(), 
                                    End = numeric())
    }
    return(error_regions)
}

check_if_blacklisted <- function(chrom_num, start_pos, end_pos, error_regions) {
    if (nrow(error_regions) == 0) {
        return(TRUE)
    }
    if (any(error_regions$Chrom == chrom_num & 
            error_regions$Start == start_pos & 
            error_regions$End == end_pos)) {
        message("Skipping blacklisted region: ", chrom_num, ":", 
                start_pos, "-", end_pos)
        return(FALSE)
    } else {
        return(TRUE) # Region is not blacklisted
    }
}

get_vmr_list <- function(region) {
    base_dir <- here("heritability/gcta", tolower(region), "_m")
    vmr_file <- here(base_dir, "vmr_list.txt")
    if (!file.exists(vmr_file)) {
        stop("VMR list file not found: ", vmr_file)
    }
    return(read.table(vmr_file, header=FALSE, stringsAsFactors=FALSE))
}

construct_data_path <- function(chrom_num, spos, epos, region, data_type) {
    chrom_dir <- paste0("chr_", chrom_num)        
    base_dir  <- here("heritability/gcta", tolower(region), "_m")
    
    if (tolower(data_type) == "plink") {
        inpath  <- "plink_format"
        data_fn <- paste0("subset_TOPMed_LIBD.AA.", spos, "_", epos, ".bed")
    } else if (tolower(data_type) == "vmr") {
        inpath  <- "vmr"
        data_fn <- paste0(spos, "_", epos, "_meth.phen")
    } else {
        stop("Unknown data_type specificed: ", data_type)
    }
    data_dir  <- here(base_dir, inpath, chrom_dir)
    data_path     <- file.path(data_dir, data_fn)
    
    if (!file.exists(data_path)) {
        stop(paste(data_type, "file not found:", data_path))
    }
    return(data_path)
}

load_genotypes <- function(geno_bed_path) {
    cat("Processing PLINK file:", basename(geno_bed_path), "\n")
                                        # Use tempfile for backingfile to
                                        # avoid conflicts in array jobs
    backing_rds <- tempfile(fileext = ".rds")
    rds_path    <- snp_readBed(geno_bed_path,
                               backingfile=sub("\\.rds$", "", backing_rds))
    return(snp_attach(rds_path))
}

load_phenotypes <- function(vmr_data_path) {
    cat("Processing VMR file:", basename(vmr_data_path), "\n")
    pheno    <- read.table(vmr_data_path, header=FALSE)
    return(pheno[, 3])
}

perform_snp_clumping <- function(G_imputed, info, pheno_scaled) {
    corrs <- big_univLinReg(G_imputed, pheno_scaled)
    stat  <- abs(corrs$estim)
    ## Default parameters are 0.2 r2 threshold and windo of 500 kb
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
region  <- Sys.getenv("region")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (is.na(task_id)) {
    stop("SLURM_ARRAY_TASK_ID is not set or is not a valid integer.")
}
if (region == "") {
    stop("Region environment variable is not set.")
}

                                        # Set reproducible seed per task
RNGkind("L'Ecuyer-CMRG")
set.seed(20250525 + task_id)

                                        # Load error regions
error_regions <- get_error_list()

                                        # Load VMR data
vmr_list      <- get_vmr_list(region)
if (task_id < 1 || task_id > nrow(vmr_list)) {
    stop("SLURM_ARRAY_TASK_ID is out of bounds for the VMR list.")
}
vmr_entry  <- vmr_list[task_id, ]
chrom_num  <- vmr_entry[[1]]
start_pos  <- vmr_entry[[2]]
end_pos    <- vmr_entry[[3]]

                                        # Check if the window is blacklisted
if (!check_if_blacklisted(chrom_num, start_pos, end_pos, error_regions)) {
    cat("Task ID ", task_id,
        " corresponds to a blacklisted region. Exiting.\n")
    quit(save = "no", status = 0) # Exit script cleanly
}

                                        # Load genotypes
geno_path <- construct_data_path(chrom_num, start_pos, end_pos, region,
                                 "PLINK")
bigSNP    <- load_genotypes(geno_path)
G         <- bigSNP$genotypes
infos     <- bigSNP$map

                                        # Load Phenotype
pheno_path <- construct_data_path(chrom_num, start_pos, end_pos, region,
                                 "VMR")
pheno_raw  <- load_phenotypes(pheno_path)
pheno_scaled <- scale(pheno_raw)[, 1] # Center and scale data

                                        # Filter out zero-variance SNPs
snp_variances <- big_apply(
    G,
    function(X, ind) {
        apply(X[, ind, drop = FALSE], 2, function(x) var(x, na.rm=TRUE))
    },
    a.combine = "c"
)
keep_idx   <- which(snp_variances > 1e-6)
infos_filt <- infos[keep_idx, ]
G_temp     <- G[, keep_idx]
G_filtered <- FBM.code256(
    nrow = nrow(G_temp), ncol = ncol(G_temp),
    code = bigSNP$genotypes$code256,
    backingfile = tempfile()
)
G_filtered[] <- G_temp[]

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
                                        # Top correlated SNPs
    if (length(residuals) != nrow(G_clumped)){
        stop("Residuals length does not match genotype matrix rows.")
    }
    
    batch_corrs <- big_univLinReg(G_clumped, residuals)
    selected_snps <- order(abs(batch_corrs$estim), decreasing = TRUE)[1:batch_size]

                                        # Fit batch via elastic net
    X_batch <- as_FBM(G_clumped[, selected_snps])
    cv_fit  <- big_spLinReg(X_batch, residuals, alphas = seq(0.05, 1, 0.05),
                            K = 5)

                                        # Extract kept indicies
    kept_ind <- attr(cv_fit, "ind.col")

    if (length(kept_ind) > 0) {
        batch_pred <- predict(cv_fit, X_batch)
        residuals  <- residuals - batch_pred
        best_betas <- summary(cv_fit, best.only = TRUE)$beta[[1]]
        global_idx <- selected_snps[kept_ind]
        
        for(i in seq_along(global_idx)){
            idx <- global_idx[i]
            accumulated_betas[1, idx] <- accumulated_betas[1, idx] + best_betas[i]
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
    if (length(lambda_idx) == 1 && !is.null(final_model$fit.preval)) {
        pred_cv   <- final_model$fit.preval[, lambda_idx]
        valid_idx <- !is.na(pheno_scaled) & !is.na(pred_cv)
        if (sum(valid_idx) > 1) {
            r_squared_cv <- cor(pheno_scaled[valid_idx], pred_cv[valid_idx])^2
        } else {
            cat("Warning: Not enough valid data points to calculate CV R^2.\n")
        }
    } else {
        cat("Warning: Using fallback prediction method for CV R^2.\n")
        pred_fallback <- predict(final_model, G_clumped[], s = "lambda.min")
        valid_idx     <- !is.na(pheno_scaled) & !is.na(pred_fallback)
        if (sum(valid_idx) > 1) {
            r_squared_cv <- cor(pheno_scaled[valid_idx], pred_fallback[valid_idx])^2
        }
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
    task_id = task_id,
    chrom = chrom_num,
    start = start_pos,
    end = end_pos,
    num_snps = ifelse(exists("G_clumped"), ncol(G_clumped), 0),
    boosting_iterations_performed = ifelse(exists("iter"), iter, 0),
    h2_unscaled = ifelse(exists("h2_unscaled"), h2_unscaled, NA),
    r_squared_cv = ifelse(exists("r_squared_cv"), r_squared_cv, NA)
)

output_df <- data.frame(
    task_id = task_id,
    chrom   = chrom_num,
    start   = start_pos,
    end     = end_pos,
    iteration = seq_along(h2_estimates),
    h2_incremental = h2_estimates
)

betas_df <- data.frame(
    task_id = task_id,
    chrom   = chrom_num,
    start   = start_pos,
    end     = end_pos,
    snp_id  = infos_filt$marker.ID,
    beta    = final_accumulated_betas
)

write.table(task_summary_df,
            file = sprintf("task_summary_stats_%d.tsv", task_id),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(output_df,
            file = sprintf("h2_estimates_%d.tsv", task_id),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(betas_df,
            file = sprintf("betas_%d.tsv", task_id),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Total SNP-based h2 (unscaled): %.4f\n", h2_unscaled))
cat(sprintf("Final h2: %.4f\n", sum(h2_estimates)))

## Reproducibility
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
