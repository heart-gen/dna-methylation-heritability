#!/usr/bin/env Rscript
library(here)
library(bigsnpr)
library(bigstatsr)

get_genotypes <- function(task_id) {
    base_dir <- here("heritability/gcta/caudate/_m")
    vmr_file <- here(base_dir, "vmr_list.txt")
                                        # Read the VMR list
    vmr_list <- read.table(vmr_file, header=FALSE, stringsAsFactors=FALSE)

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

                                        # Construct the chromosome
                                        # directory path
    chrom_dir <- paste0("chr_", chrom_num)
    geno_dir  <- here("heritability/gcta/caudate/_m/plink_format",
                      chrom_dir)

                                        # Identify corresponding PLINK file
                                        # to the VMR entry
    geno_bed_filename <- paste0("subset_TOPMed_LIBD.AA.", start_pos, "_",
                                end_pos, ".bed")
    geno_bed_path     <- file.path(geno_dir, geno_bed_filename)

                                        # Check if the PLINK file exists
    if (!file.exists(geno_bed_path)) {
        stop(paste("PLINK file not found:", geno_bed_path))
    }
    
                                        # Process the PLINK file
    cat("Processing PLINK file:", geno_bed_path, "\n")

                                        # Convert PLINK files to FBM
    geno   <- snp_readBed(geno_bed_path, backingfile = "genotypes")
    bigSNP <- snp_attach(geno)
    return(bigSNP)
}

get_vmr_data <- function(task_id) {
    # Read the VMR list
    vmr_list <- read.table("vmr_list.txt", header = FALSE,
                           stringsAsFactors = FALSE)

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

                                        # Construct the chromosome
                                        # directory path
    chrom_dir <- paste0("chr_", chrom_num)
    vmr_dir  <- here("heritability/gcta/caudate/_m/vmr", chrom_dir)

                                        # Identify corresponding VMR file
                                        # to the VMR entry
    vmr_filename <- paste0(start_pos, "_", end_pos, "_meth.phen")
    vmr_path     <- file.path(vmr_dir, vmr_filename)

                                        # Check if the VMR file exists
    if (!file.exists(vmr_path)) {
        stop(paste("VMR file not found:", vmr_path))
    }
    
                                        # Process the VMR file
    cat("Processing VMR file:", vmr_path, "\n")
    pheno <- read.table(vmr_path, header=FALSE)
    return(pheno[, 3])
}

calcalute_snp_clumping <- function(bigSNP) {
                                        # MAP-based ranking
    ind.keep <- snp_clumping(
        G = bigSNP$genotypes,
        infos.chr = bigSNP$map$chromosome,
        infos.pos = bigSNP$map$physical.pos,
        ncores = nb_cores()
    )
    return(ind.keep)
}

cal_snp_clumping_unireg <- function(bigSNP, pheno) {
                                        # Using simple univariate regression
    corrs <- big_univLinReg(bigSNP$genotypes, pheno)
    stat  <- abs(corrs$estim) # could use corrs$score, corrs$p.value, etc.
    ## Default parameters are 0.2 r2 threshold and windo of 500 kb
    ind.keep <- snp_clumping(
        G = bigSNP$genotypes,
        infos.chr = bigSNP$map$chromosome,
        infos.pos = bigSNP$map$physical.pos,
        S = stat,
        ncores = nb_cores()
    )
    return(ind.keep)
}

#### MAIN
                                        # Retrieve the SLURM array task ID
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
                                        # Load data
bigSNP  <- get_genotypes(task_id)
G       <- bigSNP$genotypes  # Returns an FBM
pheno   <- get_vmr_data(task_id)

                                        # Run clumping
clumped   <- cal_snp_clumping_unireg(bigSNP, pheno)
G_clumped <- big_copy(G, ind.col = which(clumped))

## Boosting framework
                                        # Parameters
n_iter     <- 100  # Total boosting iterations
batch_size <- 1000 # SNPs per batch
mstop      <- 10   # Batch-specific iterations

                                        # Initialize
residuals    <- pheno
h2_estimates <- numeric(n_iter)
betas        <- big_copy(G_clumped, type = "double") # Store SNP effects

                                        # Boosting loop
for (iter in 1:n_iter) {
                                        # Batch selection: Top correlated SNPs
    corrs    <- big_univLinReg(G_clumped, residuals, ncores=nb_cores())
    selected <- order(-abs(corrs$estim))[1:min(batch_size, ncol(G_clumped))]

                                        # Fit batch via elastic net
    X_batch <- big_copy(G_clumped, ind.col = selected)
    cv_fit  <- big_spLinReg(X_batch, residuals, K=5, ncores=nb_cores())

                                        # Update residuals & store effects
    pred            <- predict(cv_fit, X_batch)
    residuals       <- residuals - pred
    betas[selected] <- betas[selected] + cv_fit$beta

                                        # Calculate incremental h2
    h2_estimates[iter] <- var(pred) / var(pheno)

                                        # Early stopping
    if (iter > 10 && sd(tail(h2_estimates, 5)) < 0.001) break
}
## Heritability estimates
                                        # Final h2 estimate
final_h2 <- sum(h2_estimates, na.rm = TRUE)

                                        # SNP-specific heritability contributions
snp_h2 <- betas[]^2 * big_apply(G_clumped, function(X) apply(X, 2, var),
                                ncores = nb_cores())
total_h2 <- sum(snp_h2) / var(pheno)
