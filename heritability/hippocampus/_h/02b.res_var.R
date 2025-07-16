#### Calculate variance of residuals after regressing out PCs ####

suppressPackageStartupMessages({
  library('data.table')
  library(here)
  library(readr)
  library(stringr)
})

                                        # get chr from command line arg
args <- commandArgs(trailingOnly = TRUE)
chr  <- args[1]

# Function
remove_ct_snps <- function(f_snp, meth_levels, pos) {
  snp <- fread(f_snp, header=FALSE, data.table=FALSE)
  idx <- is.element(pos[, 1], snp[, 1])
  meth_levels <- meth_levels[, !idx]
  pos <- pos[!idx, 1]
  
  return(list(meth_levels = meth_levels, pos = pos))
}

filter_pheno <- function(meth_levels, brain_id, ances, demo, pc) {
  # keep AA only
  id2 <- intersect(intersect(ances$id[ances$group == "AA"], brain_id), 
                   demo$brnum[demo$region == "HIPPO" & demo$agedeath >= 17])
  meth_levels <- meth_levels[match(id2, brain_id), ]
  
  # align samples
  ind <- intersect(id2, pc$X)
  pc <- pc[match(ind, pc$X), -1]
  meth_levels <- meth_levels[match(ind, id2), ]
  
  return(list(meth_levels = meth_levels, pc = pc, ind = ind))
}

regress_pcs <- function(meth_levels, pc){
  pc_res <- meth_levels
  for(i in 1:ncol(meth_levels)){
    if(! i %% 100){
      cat(i,"\n")
    }
    d <- as.data.frame(cbind(y=meth_levels[,i],pc[,1:5]))
    model = lm(y ~ .,data=d)
    pc_res[,i] = resid(model)
  }
  
  return(pc_res)
}

residual_variance <- function(pc_res, pos, chr, output_path, chunk_id) {
  # get variance for residual
  res_var <- apply(pc_res, 2, sd)
  
  # write to file
  f_out <- file.path(output_path, paste0("res_var_", chunk_id, ".tsv"))
  out <- data.frame(chr=chr, pos=pos, sd=res_var)
  write.table(out, f_out, col.names=FALSE, row.names=FALSE, sep="\t", quote=F)
  return(out)
}

# Main
output_path <- here("heritability", "hippocampus", "_m", "pca", paste0("chr_", chr))
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# set file paths
cpg_meth <- here("heritability", "hippocampus", "_m", "cpg", 
                 paste0("chr_", chr), "cpg_meth.phen")
cpg_names <- here("heritability", "hippocampus", "_m", "cpg", 
                  paste0("chr_", chr), "cpg_pos.txt")
pheno_file_path <- here("inputs", "phenotypes", "_m", "phenotypes-AA.tsv")
f_ances <- here("inputs", "genetic-ancestry", 
                "structure.out_ancestry_proportion_raceDemo_compare")
f_snp <- paste0("/projects/b1213/resources/libd_data/wgbs/DEM2/snps_CT/chr", chr)
pca <- here("heritability", "hippocampus", "_m", "pca", 
            paste0("chr_", chr), "pc.csv")

# read data
pc <- read.csv(pca)
pos <- read.table(cpg_names, header=FALSE)
pos <- pos[-c(1, 2), , drop = FALSE]
ances <- read.table(f_ances, header=TRUE)
demo  <- fread(pheno_file_path, header = TRUE)

tmp_dir <- here("heritability", "hippocampus", "_m", "cpg", paste0("chr_", chr), "tmp_files")
tmp_files <- list.files(tmp_dir, pattern = "^cpg_meth_.*\\.tsv$", full.names = TRUE)
cat("Found", length(tmp_files), "files to read\n")

for (chunk_path in tmp_files) {
  cat("Processing:", chunk_path, "\n")
  flush.console()
  
  chunk_data <- fread(chunk_path, header = TRUE, data.table = FALSE)
  brain_id <- chunk_data[, 1]
  meth_levels <- chunk_data[, -c(1, 2)]
  
  # Extract CpG indices from filename
  idx <- as.integer(str_extract_all(basename(chunk_path), "\\d+")[[1]])
  start_idx <- idx[1] - 2  # Adjust for skipped rows in pos
  end_idx <- idx[2] - 2
  
  pos_chunk <- pos[start_idx:end_idx, , drop = FALSE]
  
  if (!all(colnames(meth_levels) == pos_chunk[, 1])) {
     stop("CpG names in methylation matrix do not match position file.")
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
  pc_res <- regress_pcs(pheno$meth_levels, pheno$pc)
  
  # Get chunk identifier for output
  chunk_id <- str_extract(basename(chunk_path), "(?<=cpg_meth_)[0-9]+_[0-9]+")
  
  # Write variance output
  residual_variance(pc_res, filtered$pos, chr, output_path, chunk_id)
}

res_files <- list.files(output_path, pattern = "^res_var_.*\\.tsv$", full.names = TRUE)
res_combined <- rbindlist(lapply(res_files, fread))
write.table(res_combined, file.path(output_path, "res_var_all.tsv"),
            col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

cat("Finished processing all files\n")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
