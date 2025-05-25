#!/usr/bin/env Rscript
library(here)
library(glmnet)
library(bigsnpr)
library(bigstatsr)

get_error_list <- function() {
                                        # Read error regions
    error_file <- "../_h/snp-error-window.tsv"
    if(file.exists(error_file)) {
        error_regions <- read.table(error_file, header=TRUE, sep="\t") |>
            dplyr::mutate(Chrom=as.numeric(gsub("chr", "", Chr)),
                          Start=as.numeric(Start), End=as.numeric(End)) |>
            dplyr::select(Chrom, Start, End)
    } else {
        error_regions <- data.frame(Chrom = numeric(), Start = numeric(), 
                                    End = numeric())
    }
    return(error_regions)
}

check_blacklist_regions <- function(chrom_num, start_pos, end_pos) {
                                        # Read error regions
    error_regions <- get_error_list()

                                        # Check if region is in error list
    if(any(error_regions$Chrom == chrom_num & 
           error_regions$Start == start_pos & 
           error_regions$End == end_pos)) {
        message("Skipping blacklisted region: ", chrom_num, ":", 
                start_pos, "-", end_pos)
        return(FALSE)
    } else {
        return(TRUE)
    }
}

get_vmrs <- function(region) {
    base_dir <- here("heritability/gcta", tolower(region), "_m")
    vmr_file <- here(base_dir, "vmr_list.txt")
    return(read.table(vmr_file, header=FALSE, stringsAsFactors=FALSE))
}

get_data_path <- function(chrom_num, spos, epos, region, DATA) {
    chrom_dir <- paste0("chr_", chrom_num)        
    base_dir  <- here("heritability/gcta", tolower(region), "_m")
                                        # Construct the chromosome
                                        # directory path
    if(tolower(DATA) == "plink") {
        inpath  <- "plink_format"
        data_fn <- paste0("subset_TOPMed_LIBD.AA.", spos, "_",
                          epos, ".bed")
    } else {
        inpath  <- "vmr"
        data_fn <- paste0(spos, "_", epos, "_meth.phen")
    }
    data_dir  <- here(base_dir, inpath, chrom_dir)

                                        # Identify corresponding data file
                                        # to the VMR entry
    data_path     <- file.path(data_dir, data_fn)

                                        # Check if the data file exists
    if(!file.exists(data_path)) {
        stop(paste(DATA, "file not found:", data_path))
    }
    return(data_path)
}

get_genotypes <- function(task_id, region) {
                                        # Read the VMR list
    vmr_list <- get_vmrs(region)

                                        # Validate the task ID
    if(task_id < 1 || task_id > nrow(vmr_list)) {
        stop("SLURM_ARRAY_TASK_ID is out of bounds.")
    }
                                        # Extract the specific VMR entry
                                        # for this task
    vmr_entry <- vmr_list[task_id, ]
    chrom_num <- vmr_entry[1]
    start_pos <- vmr_entry[2]
    end_pos   <- vmr_entry[3]

                                        # Check if region is in error list
    if(!check_blacklist_regions(chrom_num, start_pos, end_pos)) {
        message("Skipping task_id: ", task_id)
        return(NULL)
    }

                                        # Process the PLINK file
    geno_bed_path <- get_data_path(chrom_num, start_pos, end_pos,
                                   region, "PLINK")
    cat("Processing PLINK file:", basename(geno_bed_path), "\n")

                                        # Convert PLINK files to FBM
    geno   <- snp_readBed(geno_bed_path, backingfile = "genotypes")
    bigSNP <- snp_attach(geno)
    return(bigSNP)
}

get_vmr_data <- function(task_id, region) {
                                        # Read the VMR list
    vmr_list <- get_vmrs(region)

                                        # Validate the task ID
    if (task_id < 1 || task_id > nrow(vmr_list)) {
        stop("SLURM_ARRAY_TASK_ID is out of bounds.")
    }
                                        # Extract the specific VMR entry
                                        # for this task
    vmr_entry <- vmr_list[task_id, ]
    chrom_num <- vmr_entry[1]
    start_pos <- vmr_entry[2]
    end_pos   <- vmr_entry[3]

                                        # Check if region is in error list
    if(!check_blacklist_regions(chrom_num, start_pos, end_pos)) {
        message("Skipping task_id: ", task_id)
        return(NULL)
    }
    
                                        # Process the VMR file
    vmr_path <- get_data_path(chrom_num, start_pos, end_pos, region, "VMR")
    cat("Processing VMR file:", basename(vmr_path), "\n")
    pheno    <- read.table(vmr_path, header=FALSE)
    return(pheno[, 3])
}

cal_snp_clumping_unireg <- function(G_imputed, mapped_df, pheno) {
                                        # Using simple univariate regression
    corrs <- big_univLinReg(G_imputed, pheno)
    stat  <- abs(corrs$estim)
    ## Default parameters are 0.2 r2 threshold and windo of 500 kb
    ind.keep <- snp_clumping(
        G = G_imputed,
        infos.chr = mapped_df$chromosome,
        infos.pos = mapped_df$physical.pos,
        S = stat
    )
    return(ind.keep)
}

#### MAIN
                                        # Retrieve variables
region  <- Sys.getenv("region")
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

                                        # Set reproducible seed per task
RNGkind("L'Ecuyer-CMRG")
set.seed(20250525 + task_id)

                                        # Load VMR info
vmr_list   <- get_vmrs(region)
vmr_entry  <- vmr_list[task_id, ]
chrom_num  <- vmr_entry[1]
start_pos  <- vmr_entry[2]
end_pos    <- vmr_entry[3]

                                        # Load data
bigSNP  <- get_genotypes(task_id, region)
pheno   <- get_vmr_data(task_id, region)
pheno   <- scale(pheno)[, 1] # center and scale
infos   <- bigSNP$map

# Filter out zero-variance SNPs
variances <- big_apply(
    bigSNP$genotypes,
    function(X, ind) {
        apply(X[, ind, drop = FALSE], 2, function(x) var(x, na.rm=TRUE))
    },
    a.combine = "c"
)
keep       <- variances > 0
infos      <- infos[keep, ]

G_filtered <- as_FBM(bigSNP$genotypes[, keep])

                                        # Impute missing values
G_imputed <- snp_fastImputeSimple(G_filtered, method = "mode")

                                        # Run clumping
clumped   <- cal_snp_clumping_unireg(G_imputed, infos, pheno)
G_clumped <- as_FBM(G_imputed[, clumped])

## Boosting framework
                                        # Boosting parameters
n_iter     <- 100  # Total boosting iterations
batch_size <- 1000 # SNPs per batch
mstop      <- 10   # Batch-specific iterations

                                        # Initialize
residuals    <- pheno
h2_estimates <- numeric(n_iter)
betas        <- big_copy(G_clumped, type = "double") # Store SNP effects

# Boosting loop
for (iter in 1:n_iter) {
                                        # Batch selection:
                                        # Top correlated SNPs
    corrs    <- big_univLinReg(G_clumped, residuals)
    selected <- order(-abs(corrs$estim))[1:min(batch_size, ncol(G_clumped))]

                                        # Fit batch via elastic net
    X_batch <- as_FBM(G_clumped[, selected])
    cv_fit  <- big_spLinReg(X_batch, residuals, alphas = seq(0.05, 1, 0.05),
                            K = 5)

                                        # Extract kept indicies
    kept_ind <- attr(cv_fit, "ind.col")

                                        # Update residuals & store effects
    pred            <- predict(cv_fit, X_batch)
    residuals       <- residuals - pred

    best_betas      <- summary(cv_fit, best.only = TRUE)$beta[[1]]
    global_kept_ind <- selected[kept_ind]
    betas[global_kept_ind] <- betas[global_kept_ind] + best_betas

                                        # Calculate incremental h2
    h2_estimates[iter] <- var(pred)

                                        # Early stopping
    if (iter > 10 && sd(tail(h2_estimates, 5)) < 0.001) break
}

## Refit using OLS or Ridge and evaluate R^2
G_scaled    <- big_scale()(G_clumped)$scale

final_model <- cv.glmnet(G_scaled[], pheno, alpha = 0, nfolds = 5)
pred_cv     <- predict(final_model, G_scaled[], s = "lambda.min")
r2_cv       <- cor(pheno, pred_cv)^2

## Heritability estimates (standardized)
var_y    <- 1  # because pheno was scaled
h2_final <- sum(betas[]^2)

                                        # SNP-wise vars, unscaled estimate
snp_vars <- big_apply(
    G_clumped, function(X, ind) {
        apply(X[, ind, drop = FALSE], 2, var)
    },
    a.combine = "c"
)
h2_unscaled <- sum(betas[]^2 * snp_vars) / var(pheno)

                                        # Save h2_estimates and betas
output_df <- data.frame(
    task_id = task_id,
    chrom   = chrom_num,
    start   = start_pos,
    end     = end_pos,
    iteration = seq_along(h2_estimates),
    h2_incremental = h2_estimates
)

betas_df <- data.frame(
    snp_index = seq_along(betas),
    beta      = betas[],
    chrom     = chrom_num,
    start     = start_pos,
    end       = end_pos,
    task_id   = task_id
)

write.table(output_df,
            file = sprintf("h2_estimates_task_%d.tsv", task_id),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(betas_df,
            file = sprintf("betas_task_%d.tsv", task_id),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Final cumulative h2 (standardized): %.4f\n", h2_final))
cat(sprintf("Total SNP-based h2 (unscaled): %.4f\n", h2_unscaled))
cat(sprintf("Cross-validated R^2: %.4f\n", r2_cv))

## Reproducibility
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
