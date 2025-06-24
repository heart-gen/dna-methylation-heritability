#### Perform principal component analysis ####

suppressPackageStartupMessages({
    library('bsseq')
    library('HDF5Array')
    library(DelayedMatrixStats)
    library('data.table')
    library(here)
})

                                        # get chr from command line arg
args <- commandArgs(trailingOnly = TRUE)
chr  <- args[1]

## Function
get_top_meth <- function(BSobj, v) {

                                        # get top 1M variable CpG
    v_top  <- v[order(v$sd * -1), ]
    v_top  <- v_top[1:10^6, ]

                                        # DNAm levels for top 1M
    tmp    <- v_top[v_top$chr == chr, ]
    BS     <- BSobj[is.element(start(BSobj), tmp$start), ]
    meth   <- as.matrix(getMeth(BS))
    id_top <- colData(BS)$brnum
    colnames(meth) <- id_top
    return(list(BS = BS, meth = meth, id_top = id_top))
}

get_ances <- function(ances_file_path, id_top) {
    gen_ances <- read.table(ances_file_path, header=TRUE)
    idx       <- match(id_top, gen_ances$id)
    return(list(gen_ances = gen_ances, idx = idx))
}

res_ances <- function(meth, gen_ances, idx) {
    res <- meth
    for(i in 1:nrow(meth)) {
      if (i %% 10000 == 0) {
        cat(i, "\n")
      }
      model    <- lm(meth[i, ] ~ gen_ances$Afr[idx])
      res[i, ] <- resid(model)
    }
    return(res)
}

pca <- function(res, output_path) {
    pc <- prcomp(t(res), scale. = TRUE, center = TRUE)
    print(cumsum(summary(pc)$importance[2,1:10]))

                                        # plot PCA
    pdf(file = file.path(output_path, "pca.pdf"))
    plot(pc)
    plot(pc$x[, 1], pc$x[, 2])
    dev.off()

                                        # write file
    write.csv(pc$x, file = file.path(output_path, "pc.csv"))
    return(pc)
}

corr_pc_ances <- function(pc, gen_ances, idx, output_path) {
    res <- matrix(NA, nrow=ncol(pc$x), ncol=2)
    colnames(res) <- c("Correlation", "p-value")
    for(i in 1:ncol(pc$x)){
      tmp       <- cor.test(gen_ances$Afr[idx], pc$x[, i])
      res[i, 1] <- tmp$estimate
      res[i, 2] <- tmp$p.value
    }
    write.csv(res, file = file.path(output_path, "pc_ances_cor.csv"))
}

remove_ct_snps <- function(snp_file_path, meth, pos) {
    snp <- fread(snp_file_path, header=F, data.table=F)
    idx <- is.element(pos[,1], snp[,1])
    meth <- meth[,!idx]
    pos <- pos[!idx,1]
    return(list(meth, pos))
}

filter_pheno <- function(brnum, meth, pc, ances_file_path, pheno_file_path) {
  ances <- read.table(ances_file_path, header=TRUE)
  demo  <- read.csv(pheno_file_path, header = TRUE)
  samples <- intersect(intersect(ances$brnum[ances$group == "AA"],brnum),
                       demo$BrNum[demo$Region == "Caudate" & demo$Age >= 17])
  meth <- meth[match(samples,brnum),]
  
  # align samples
  ind <- intersect(samples, pc$X)
  pc <- pc[match(ind, pc$X), -1]
  meth <- meth[match(ind, samples), ]
  return(list(meth, pc))
}

regress_pcs <- function(meth, pc, chr, pos, output_path){
  residuals <- meth
  for(i in 1:ncol(meth)){
    if(! i %% 100){
      cat(i,"\n")
    }
    d <- as.data.frame(cbind(y=meth[,i],pc[,1:5]))
    model = lm(y ~ .,data=d)
    residuals[,i] = resid(model)
  }
  
  # get variance for residual
  res_var <- apply(residuals,2,sd)
  
  # output
  f_out <- file.path(output_path, "res_var.tsv")
  out <- data.frame(chr=chr, pos=pos, sd=res_var)
  write.table(out, f_out, col.names=F, row.names=F, sep="\t", quote=F)
  return(out)
}


## Main
                                        # create output directories if they  
                                        # don't exist
output_path <- here("heritability", "gcta", "caudate", "_m", "pca", paste0("chr_", chr))
if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

                                        # load sd of raw DNAm
load(here("heritability", "gcta", "caudate", "_m", "cpg", paste0("chr_", chr), "stats.rda"))
v <- data.frame(chr = chr, start = start(BSobj), sd = sds)

# set file paths
cpg_meth <- here("heritability", "gcta", "caudate", "_m", "cpg", paste0("chr_", chr), "cpg_meth.phen")
cpg_names <- here("heritability", "gcta", "caudate", "_m", "cpg", paste0("chr_", chr), "cpg_pos.txt")
pheno_file_path <- here("inputs/phenotypes/merged/_m/merged_phenotypes.csv")
snp_file_path <- paste0("/projects/b1213/resources/libd_data/wgbs/DEM2/snps_CT/chr",chr)

# read data
meth <- fread(cpg_meth, header=TRUE, data.table=FALSE)
brnum <- meth[, 1]
meth <- meth[,-c(1, 2)]
pos <- read.table(cpg_names, header=FALSE)
pos <- pos[-c(1, 2), , drop = FALSE]

meth_top <- get_top_meth(BSobj, v)

                                        # load genetic ancestry
ances_file_path <- here("inputs/genetic-ancestry/structure.out_ancestry_proportion_raceDemo_compare")
ances   <- get_ances(ances_file_path, meth_top$id_top)

                                        # residual after regressing out
                                        # Afr proportion
res <- res_ances(meth_top$meth, ances$gen_ances, ances$idx)

                                        # perform and plot pca
pc <- pca(res, output_path)
#PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10 
#0.05648 0.08545 0.11074 0.13106 0.15105 0.16687 0.17921 0.19071 0.20140 0.21163 

                                        # correlation of pc with ancestry
corr_pc_ances(pc, ances$gen_ances, ances$idx, output_path)

filtered_cpg <- remove_ct_snps(snp_file_path, meth, pos)

filtered_samples <- filter_pheno(brnum, filtered_cpg$meth, pc, ances_file_path, pheno_file_path)

out <- regress_pcs(filtered_samples$meth, filtered_samples$pc, chr, filtered_cpg$pos, output_path)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
