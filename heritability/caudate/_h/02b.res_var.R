#### Calculate variance of residuals after regressing out PCs ####

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(stringr)
  library(data.table)
  library(matrixStats)
})

                                        # Get chr from command line arg
args <- commandArgs(trailingOnly = TRUE)
chr  <- args[1]

## Function
remove_ct_snps <- function(f_snp, meth_levels, pos) {
    snp <- fread(f_snp, header = FALSE, data.table = FALSE)[, 1]
    keep_idx <- !pos[[1]] %in% snp
    meth_levels <- meth_levels[, keep_idx, drop = FALSE]
    pos <- pos[keep_idx, , drop = FALSE]
    return(list(meth_levels = meth_levels, pos = pos))
}

filter_pheno <- function(meth_levels, brain_id, ances, demo, pc) {
                                        # Filter ancestry and phenotypes
    demo  <- demo[region == "caudate" & agedeath >= 17]
    ances <- ances[group == "AA"]

                                        # Keep AA only
    valid_ids <- intersect(intersect(ances$id, demo$brnum), brain_id)
    valid_ids <- intersect(valid_ids, pc$V1)

                                        # Align samples
    meth_levels <- meth_levels[match(valid_ids, brain_id), , drop = FALSE]
    pc_filt     <- pc[V1 %in% valid_ids, -1, with = FALSE]
    return(list(meth_levels = meth_levels, pc = pc_filt, ind = valid_ids))
}

## regress_pcs <- function(meth_levels, pc){
##     pc_res  <- meth_levels
##     for(i in 1:ncol(meth_levels)){
##         if(! i %% 100){
##             cat(i,"\n")
##         }
##         temp_df    <- data.frame(y = meth_levels[, i], pc[, 1:5]))
##         model      <- lm(y ~ .,data=temp_df)
##         pc_res[,i] <- resid(model)
##     }
##     return(pc_res)
## }

regress_pcs_vectorized <- function(meth_levels, pc_matrix) {
    meth_t <- t(meth_levels)
    design <- cbind(1, as.matrix(pc_matrix[, 1:5, drop = FALSE]))
    fit    <- limma::lmFit(meth_t, design)
    resids <- limma::residuals.MArrayLM(fit, meth_t)
    return(t(resids))
}

residual_variance <- function(pc_res, pos, chr, output_path, chunk_id) {
                                        # Get variance for residual
    res_var <- colVars(pc_res) # Add variance
    res_sd  <- colSds(pc_res)
    out     <- data.frame(chr = chr, pos = pos[[1]], sd = res_sd, var = res_var)
    return(out)
}

## Main
output_path <- here("heritability", "caudate", "_m", "pca", paste0("chr_", chr))
if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

                                        # Set file paths
cpg_names <- here("heritability", "caudate", "_m", "cpg",
                  paste0("chr_", chr), "cpg_pos.txt")
pheno_file_path <- here("inputs", "phenotypes", "_m", "phenotypes-AA.tsv")
f_ances <- here("inputs", "genetic-ancestry",
                "structure.out_ancestry_proportion_raceDemo_compare")
f_snp <- paste0("/projects/b1213/resources/libd_data/wgbs/DEM2/snps_CT/chr",chr)
pca <- here("heritability", "caudate", "_m", "pca",
            paste0("chr_", chr), "pc.csv")

                                        # Load data
pc    <- fread(pca)
pos   <- fread(cpg_names, header = FALSE)[-c(1, 2), , drop = FALSE]
ances <- fread(f_ances,)
demo  <- fread(pheno_file_path)

tmp_dir   <- here("heritability", "caudate", "_m", "cpg", paste0("chr_", chr),
                  "tmp_files")
tmp_files <- list.files(tmp_dir, pattern = "^cpg_meth_.*\\.tsv$",
                        full.names = TRUE)
cat("Found", length(tmp_files), "files to read\n")

                                        # Process each chunk
output_file <- file.path(output_path, "res_var_all.tsv")
first_chunk <- TRUE

for (chunk_path in tmp_files) {
    cat("Processing:", chunk_path, "\n"); flush.console()

    chunk_data  <- fread(chunk_path, header = TRUE)
    brain_id    <- chunk_data[[1]]
    meth_levels <- as.matrix(chunk_data[, -c(1, 2), with = FALSE])

                                        # Extract CpG indices from filename
    idx <- as.integer(str_extract_all(basename(chunk_path), "\\d+")[[1]])
    pos_chunk <- pos[(idx[1] - 2):(idx[2] - 2), , drop = FALSE]

    if (!identical(colnames(meth_levels), as.character(pos_chunk[[1]]))) {
        stop("CpG names in methylation matrix do not match position file: ",
             basename(chunk_path))
    }

    if (ncol(meth_levels) != nrow(pos_chunk)) {
        stop("Mismatch between meth_levels and pos_chunk: ",
             ncol(meth_levels), " vs ", nrow(pos_chunk))
    }

                                        # Remove CT SNPs
    filtered <- remove_ct_snps(f_snp, meth_levels, pos_chunk)

                                        # Filter and align phenotype
    pheno <- filter_pheno(filtered$meth_levels, brain_id, ances, demo, pc)

                                        # Regress PCs
    pc_res <- regress_pcs_vectorized(pheno$meth_levels, pheno$pc)

                                        # Get chunk identifier for output
    chunk_id <- str_extract(basename(chunk_path), "(?<=cpg_meth_)[0-9]+_[0-9]+")

                                        # Calculate variances
    out <- residual_variance(pc_res, filtered$pos, chr, output_path, chunk_id)
    fwrite(
        out, output_file, append =!first_chunk,
        col.names = first_chunk, sep = "\t", quote = FALSE)
    first_chunk <- FALSE
    rm(chunk_data, meth_levels, pc_res, out); gc()
}

cat("Finished processing all files\n")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
