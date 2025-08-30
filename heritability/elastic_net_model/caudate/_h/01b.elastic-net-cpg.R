## --- Load libraries --- ##
library(here)
library(glmnet)
library(bigsnpr)
library(bigstatsr)
library(data.table)

## --- Helper Functions --- ##
#get_cpg_list <- function(region) {
#    base_dir <- file.path("/projects/b1042/HEART-GeN-Lab/dna_meth_heritability/cpg_snps", 
#                          paste0(tolower(region), "/_m/cpg/chr_21/test/chunked_cpg"))
#    cpg_file <- file.path(base_dir, "chunk_000.bed")
#    if (!file.exists(cpg_file)) {
#        stop("CPG list file not found: ", cpg_file)
#    }
#    return(read.table(cpg_file, header=FALSE, stringsAsFactors=FALSE))
#}

get_sample_list <- function(region) {
    base_dir <- here("heritability", tolower(region), "_m")
    sample_file <- here(base_dir, "samples.txt")
    if (!file.exists(sample_file)) {
        stop("Sample list file not found: ", sample_file)
    }
    return(sample_file)
}

construct_data_path <- function(chrom_num, spos, epos, region, data_type) {
    chrom_dir <- paste0("chr_", chrom_num)
    base_dir  <- here("heritability", tolower(region), "_m")

    if (tolower(data_type) == "cpg") {
        inpath  <- "cpg"
        data_fn <- paste0("res_cpg_meth.phen")
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

extract_genotypes <- function(plink_base, samples, chrom, start, end) {
    temp_dir <- tempdir()
    out_prefix <- tempfile(tmpdir = temp_dir)
                                            # extract SNPs using plink

    system2("/software/plink/plink_2.0_alpha_3.3_x86_64/plink2", args = c(
        "--pfile", plink_base,
        "--chr", chrom,
        "--from-bp", start,
        "--to-bp", end,
        "--make-bed",
        "--keep", samples,
        "--no-parents",
        "--no-sex",
        "--no-pheno",
        "--out", out_prefix
    ))

    bedfile <- paste0(out_prefix, ".bed")
    if (!file.exists(bedfile)) stop("PLINK extraction failed: .bed file not found.")

    cat("Processing PLINK file:", basename(bedfile), "\n", flush = TRUE)
                                        # Use tempfile for backingfile to
                                        # avoid conflicts in array jobs
    backing_rds <- tempfile(fileext = ".rds")
    rds_path    <- snp_readBed(bedfile,
                               backingfile=sub("\\.rds$", "", backing_rds))
    file.remove(paste0(out_prefix, c(".bed", ".bim", ".fam", ".log")))
    
    return(snp_attach(rds_path))
}

match_cpgs_to_chunk <- function(chunk_regions, cpg_df, phenotype_file) {
  colnames(chunk_regions) <- c("chr", "start", "end")
  
                                        # Get CpG site position
  chunk_regions$start_adj <- chunk_regions$start + 20000
  chunk_regions$end_adj   <- chunk_regions$end - 20000
  
                                        # Match CpGs by position
  overlaps <- do.call(rbind, lapply(1:nrow(chunk_regions), function(i) {
    subset(cpg_df, pos >= chunk_regions$start_adj[i] & pos <= chunk_regions$end_adj[i])
  }))
  
  if (nrow(overlaps) == 0) {
    stop("No CpGs matched in chunk: ", chunk_bed_file)
  }
  
  cpg_cols <- sort(unique(overlaps$col_idx))
  cols_to_keep <- c(1, 2, cpg_cols)  # FID, IID + matched CpGs
  
                                          # Read subset 
  phenos_subset <- fread(phenotype_file, header = TRUE, 
                         select = cols_to_keep, data.table = FALSE)
  
  return(list(
    phenos = phenos_subset,
    matched_cpgs = overlaps
  ))
}

#load_phenotypes <- function(cpg_data_path) {
#    cat("Processing CPG file:", basename(cpg_data_path), "\n")
#    pheno    <- read.table(cpg_data_path, header=FALSE)
#    return(pheno[, 3])
#}

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
task_id <- as.integer(Sys.getenv("task_id"))

if (is.na(task_id)) {
    stop("SLURM_ARRAY_TASK_ID is not set or is not a valid integer.")
}
if (region == "") {
    stop("Region environment variable is not set.")
}

                                        # Set reproducible seed per task
RNGkind("L'Ecuyer-CMRG")
set.seed(20250525 + task_id)

                                        # Read CpG mapping
cpg_pos_file <- here("heritability", tolower(region), "_m", "cpg", "chr_21", "cpg_pos.txt")
pos <- as.integer(readLines(cpg_pos_file)[grepl("^\\d+$", readLines(cpg_pos_file))])

cpg_df <- data.frame(
  pos = pos,
  col_idx = seq_along(pos) + 2  # account for FID, IID
)

base_dir <- file.path("/projects/b1042/HEART-GeN-Lab/dna_meth_heritability/cpg_snps", 
                      paste0(tolower(region), "/_m/cpg/chr_21/test/chunked_cpg"))
chunk_bed_file <- file.path(base_dir, "chunk_1000.bed")
chunk_regions <- read.table(chunk_bed_file, header = FALSE, stringsAsFactors = FALSE)
colnames(chunk_regions) <- c("chr", "start", "end")

                                          # Load sample list
samples <- get_sample_list(region)

phenotype_file <- here("heritability", paste0(tolower(region)), "_m", 
                       "cpg", "chr_21", "res_cpg_meth.phen")

                                          # Load Phenotype
pheno_result <- match_cpgs_to_chunk(chunk_regions, cpg_df, phenotype_file)
pheno_raw <- pheno_result$phenos[, -c(1,2)]  # Drop FID/IID
matched_cpgs <- pheno_result$matched_cpgs

                                        # Load cpg data
#cpg_list      <- get_cpg_list(region)
if (task_id < 1 || task_id > nrow(chunk_regions)) {
    stop("SLURM_ARRAY_TASK_ID is out of bounds for the VMR list.")
}

num_regions <- nrow(chunk_regions)

for (i in seq_len(num_regions)) {
    chrom_num  <- chunk_regions$chr[i]
    start_pos  <- chunk_regions$start[i]
    end_pos    <- chunk_regions$end[i]
  
    pheno_scaled <- scale(pheno_raw)[, i] # Center and scale data
    
                                        # Extract genotypes
    plink_base <- here("inputs", "genotypes", "TOPMed_LIBD.AA")
    bigSNP <- extract_genotypes(plink_base, samples, chrom_num, start_pos, end_pos)
    G         <- bigSNP$genotypes
    infos     <- bigSNP$map
    
    
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
    cat("Imputing missing genotypes using mode...\n", flush = TRUE)
    G_imputed <- snp_fastImputeSimple(G_filtered, method = "mode")
    
                                        # Run SNP clumping
    clumped_idx <- perform_snp_clumping(G_imputed, infos_filt, pheno_scaled)
    G_clumped   <- as_FBM(G_imputed[, clumped_idx])
    infos_filt  <- infos_filt[clumped_idx, ]
    cat("Number of SNPs before clumping: ", ncol(G_imputed), "\n", flush = TRUE)
    cat("Number of SNPs after clumping: ", ncol(G_clumped), "\n", flush = TRUE)
    
    if (ncol(G_clumped) == 0) {
      cat("No SNPs left after clumping. Exiting.\n")
      # Write empty results or handle as needed
      quit(save = "no", status = 0)
    }
    
    ## --- Boosting framework --- ##
    cat("Starting boosting framework...\n", flush = TRUE)
    n_iter     <- 100  # Total boosting iterations
    batch_size <- min(1000, ncol(G_clumped)) # SNPs per batch
    
                                          # Initialize
    residuals    <- pheno_scaled
    h2_estimates <- numeric(n_iter)
    accumulated_betas <- FBM(1, ncol(G_clumped), type = "double", init = 0,
                             backingfile = tempfile())
    
                                          # Boosting loop
    for (iter in 1:n_iter) {
      cat("Boosting iteration: ", iter, "\n", flush = TRUE)
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
        
        for(j in seq_along(global_idx)){
          idx <- global_idx[j]
          accumulated_betas[1, idx] <- accumulated_betas[1, idx] + best_betas[j]
        }
                                          # Calculate incremental h2
        h2_estimates[iter] <- var(batch_pred)
      } else {
        h2_estimates[iter] <- 0
      }
      cat("Incremental h2 this iteration: ",
          sprintf("%.5f", h2_estimates[iter]), "\n", flush = TRUE)
      
                                          # Early stopping condition
      if (iter > 10 && sd(tail(h2_estimates[1:iter], 5)) < 0.0001) {
        cat("Early stopping criterion met at iteration: ", iter, "\n", flush = TRUE)
        h2_estimates <- h2_estimates[1:iter] # Trim unused part
        break
      }
    }
    
    ## --- Final Model Refit and Evaluation --- ##
    r_squared_cv <- NA
    if (ncol(G_clumped) > 0) {
      cat("Refitting final model using Ridge regression...\n", flush = TRUE)
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
                  r_squared_cv), flush = TRUE)
    } else {
      cat("Skipping Ridge regression as no SNPs are available.\n", flush = TRUE)
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
    
    region_idx <- i
    out_idx <- sprintf("chunk_%03d_region_%03d", task_id, region_idx)
    dir.create("cpg_summary", recursive = TRUE, showWarnings = FALSE)
    write.table(task_summary_df,
                file = file.path("cpg_summary", paste0("task_summary_stats_", out_idx, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    
    dir.create("cpg_h2", recursive = TRUE, showWarnings = FALSE)
    write.table(output_df,
                file = file.path("cpg_h2", paste0("h2_estimates_", out_idx, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    dir.create("cpg_betas", recursive = TRUE, showWarnings = FALSE)
    write.table(betas_df,
                file = file.path("cpg_betas", paste0("betas_", out_idx, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat(sprintf("Total SNP-based h2 (unscaled): %.4f\n", h2_unscaled), flush = TRUE)
    cat(sprintf("Final h2: %.4f\n", sum(h2_estimates)), flush = TRUE)
    
}

                                        # Load Phenotype
#pheno_path <- construct_data_path(chrom_num, start_pos, end_pos, region,
#                                 "cpg")
#pheno_raw  <- load_phenotypes(pheno_path)
#pheno_scaled <- scale(pheno_raw)[, 1] # Center and scale data

## Reproducibility
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()