#### Calculate variance of residuals after regressing out PCs ####

suppressPackageStartupMessages({
  library('data.table')
  library(here)
})

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
                   demo$BrNum[demo$Region == "Caudate" & demo$Age >= 17])
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

residual_variance <- function(pc_res, pos, chr, output_path) {
  # get variance for residual
  res_var <- apply(pc_res, 2, sd)
  
  # write to file
  f_out <- file.path(output_path, "res_var.tsv")
  out <- data.frame(chr=chr, pos=pos, sd=res_var)
  write.table(out, f_out, col.names=FALSE, row.names=FALSE, sep="\t", quote=F)
  return(out)
}

# Main
output_path <- here("heritability", "gcta", "caudate", "_m", "pca", paste0("chr_", chr))
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# set file paths
cpg_meth <- here("heritability", "gcta", "caudate", "_m", "cpg", 
                 paste0("chr_", chr), "cpg_meth.phen")
cpg_names <- here("heritability", "gcta", "caudate", "_m", "cpg", 
                  paste0("chr_", chr), "cpg_pos.txt")
pheno_file_path <- here("inputs", "phenotypes", "merged", 
                        "_m", "merged_phenotypes.csv")
f_ances <- here("inputs", "genetic-ancestry", 
                "structure.out_ancestry_proportion_raceDemo_compare")
f_snp <- paste0("/projects/b1213/resources/libd_data/wgbs/DEM2/snps_CT/chr", chr)
pca <- here("heritability", "gcta", "caudate", "_m", "pca", 
            paste0("chr_", chr), "pc.csv")

# read data
pc <- read.csv(pca)
meth_levels <- fread(cpg_meth, header=TRUE, data.table=FALSE)
brain_id <- meth_levels[, 1]
meth_levels <- meth_levels[, -c(1, 2)]
pos <- read.table(cpg_names, header=FALSE)
pos <- pos[-c(1, 2), , drop = FALSE]
ances <- read.table(f_ances, header=TRUE)
demo  <- read.csv(pheno_file_path, header = TRUE)

# remove CT snps at CpG sites
filtered_cpg <- remove_ct_snps(f_snp, meth_levels, pos)

# filter and align samples
filtered_samples <- filter_pheno(filtered_cpg$meth_levels, brain_id, ances, demo, pc)
  
# regress out top5 PC
pc_res <- regress_pcs(filtered_samples$meth_levels, filtered_samples$pc)

# calculate variance and write to file
res_out <- residual_variance(pc_res, filtered_cpg$pos, chr, output_path)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()