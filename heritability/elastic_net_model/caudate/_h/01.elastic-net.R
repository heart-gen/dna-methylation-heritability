#!/usr/bin/env Rscript
library(here)
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

get_data_path <- function(chrom_num, region, DATA) {
    chrom_dir <- paste0("chr_", chrom_num)        
    base_dir  <- here("heritability/gcta", tolower(region), "_m")
                                        # Construct the chromosome
                                        # directory path
    if(tolower(DATA) == "plink") {
        inpath  <- "plink_format"
        data_fn <- paste0("subset_TOPMed_LIBD.AA.", start_pos, "_",
                          end_pos, ".bed")
    } else {
        inpath  <- "vmr"
        data_fn <- paste0(start_pos, "_", end_pos, "_meth.phen")
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

                                        # Read error regions
    error_regions <- get_error_list()
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
    geno_bed_path <- get_data_path(chrom_num, region, "PLINK")
    cat("Processing PLINK file:", geno_bed_path, "\n")

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
    vmr_path <- get_data_path(chrom_num, region, "VMR")
    cat("Processing VMR file:", vmr_path, "\n")
    pheno    <- read.table(vmr_path, header=FALSE)
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
