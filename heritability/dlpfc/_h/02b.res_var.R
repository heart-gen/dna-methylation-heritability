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
    demo  <- demo[region == "dlpfc" & agedeath >= 17]
    ances <- ances[group == "AA"]

                                        # Keep AA only
    valid_ids <- intersect(intersect(ances$id, demo$brnum), brain_id)
    valid_ids <- intersect(valid_ids, pc$V1)

                                        # Align samples
    meth_levels <- meth_levels[match(valid_ids, brain_id), , drop = FALSE]
    pc_filt     <- pc[V1 %in% valid_ids, -1, with = FALSE]
    return(list(meth_levels = meth_levels, pc = pc_filt, ind = valid_ids))
}

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
output_path <- here("heritability", "dlpfc", "_m", "pca", paste0("chr_", chr))
if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

                                        # Set file paths
cpg_names <- here("heritability", "dlpfc", "_m", "cpg",
                  paste0("chr_", chr), "cpg_pos.txt")
pheno_file_path <- here("inputs", "phenotypes", "_m", "phenotypes-AA.tsv")
f_ances <- here("inputs", "genetic-ancestry",
                "structure.out_ancestry_proportion_raceDemo_compare")
f_snp <- paste0("/projects/b1213/resources/libd_data/wgbs/DEM2/snps_CT/chr",chr)
pca <- here("heritability", "dlpfc", "_m", "pca",
            paste0("chr_", chr), "pc.csv")
psam_file <- here("inputs/genotypes/TOPMed_LIBD.AA.psam")

                                        # Load data
pc        <- fread(pca)
pos       <- fread(cpg_names, header = FALSE)[-c(1, 2), , drop = FALSE]
ances     <- fread(f_ances)
demo      <- fread(pheno_file_path)
samples   <- fread(psam_file, header = FALSE, 
                   col.names = c("FID", "IID", "PAT"))[, .(FID, IID)]


tmp_dir   <- here("heritability", "dlpfc", "_m", "cpg", paste0("chr_", chr),
                  "tmp_files")
tmp_files <- list.files(tmp_dir, pattern = "^cpg_meth_.*\\.tsv$",
                        full.names = TRUE)
cat("Found", length(tmp_files), "files to read\n")

                                        # Process each chunk
output_file <- file.path(output_path, "res_var_all.tsv")
first_chunk <- TRUE

cpg_path <- here("heritability", "dlpfc", "_m", "cpg", paste0("chr_", chr))
res_meth_file <- file.path(cpg_path, paste0("res_cpg_meth.phen"))

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
    filtered_samples <- samples[match(pheno$ind, samples$FID), c("FID", "IID")]
    res_chunk <- data.table::data.table(FID = filtered_samples$FID,
                                        IID = filtered_samples$IID,
                                        pc_res)
    colnames(res_chunk) <- c("FID", "IID", as.character(filtered$pos[[1]]))
    
                                        # Get chunk identifier for output
    chunk_id <- str_extract(basename(chunk_path), "(?<=cpg_meth_)[0-9]+_[0-9]+")
    
                                        # Write intermediate residuals
    out_file_chunk <- file.path(tmp_dir, paste0("residuals_", chunk_id, ".tsv"))
    fwrite(res_chunk, out_file_chunk, sep = "\t", quote = FALSE, col.names = TRUE)
    
                                        # Calculate variances
    out <- residual_variance(pc_res, filtered$pos, chr, output_path, chunk_id)
    fwrite(
        out, output_file, append =!first_chunk,
        col.names = first_chunk, sep = "\t", quote = FALSE)
    first_chunk <- FALSE
    rm(chunk_data, meth_levels, pc_res, out); gc()
}

                                        # List all residual files
chunk_files <- list.files(tmp_dir, pattern = "^residuals_.*\\.tsv$", full.names = TRUE)
res_chunks_list <- lapply(chunk_files, function(f) fread(f, header = TRUE))
ref_samples <- res_chunks_list[[1]][, .(FID, IID)]

                                        # Read and merge by FID/IID
residuals_only <- lapply(res_chunks_list, function(dt) dt[, -c("FID", "IID"), with = FALSE])
combined_residuals <- cbind(ref_samples, do.call(cbind, residuals_only))

                                        # Write combined matrix
fwrite(combined_residuals, res_meth_file, sep = "\t", quote = FALSE)

cat("Finished processing all files\n")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
