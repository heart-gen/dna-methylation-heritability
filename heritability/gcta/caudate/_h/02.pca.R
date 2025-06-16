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

get_ances <- function(f_ances, id_top) {
    gen_ances <- read.table(f_ances, header=TRUE)
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
    return(res)
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

########
cpg_meth <- here("heritability", "gcta", "caudate", "_m", "cpg", paste0("chr_", chr), "cpg_meth.phen")
cpg_names <- here("heritability", "gcta", "caudate", "_m", "cpg", paste0("chr_", chr), "cpg_pos.txt")
f_ances <- here("inputs/genetic-ancestry/structure.out_ancestry_proportion_raceDemo_compare")
pheno_file_path <- here("inputs/phenotypes/merged/_m/merged_phenotypes.csv")
pca <- here("heritability", "gcta", "caudate", "_m", "pca", paste0("chr_", chr), "pc.csv")
pc <- read.csv(pca)

# read data - need to fix header
p <- fread(cpg_meth, header=TRUE, data.table=FALSE)
id <- p[, 1]
p <- p[,-c(1, 2)]
p_names <- read.table(cpg_names, header=FALSE)
p_names <- p_names[-c(1, 2), , drop = FALSE]
ances <- read.table(f_ances, header=TRUE)
demo  <- read.csv(pheno_file_path, header = TRUE)

# remove CT snps at CpG sites
f_snp <- paste0("/projects/b1213/resources/libd_data/wgbs/DEM2/snps_CT/chr",chr)
snp <- fread(f_snp,header=F,data.table=F)
idx <- is.element(p_names[,1], snp[,1])
p <- p[,!idx]
p_names <- p_names[!idx,1]

# keep AA only
id2 <- intersect(intersect(ances$id[ances$group == "AA"],id),demo$BrNum[demo$Region == "Caudate" & demo$Age >= 17])
p <- p[match(id2,id),]

# align samples
ind <- intersect(id2,pc$X)
pc <- pc[match(ind,pc$X),-1]
p <- p[match(ind,id2),]

# regress out top5 PC
res <- p
for(i in 1:ncol(p)){
  if(! i %% 100){
    cat(i,"\n")
  }
  d <- as.data.frame(cbind(y=p[,i],pc[,1:5]))
  model = lm(y ~ .,data=d)
  res[,i] = resid(model)
}

# get variance for residual
res_var <- apply(res,2,sd)

# output
f_out <- file.path(output_path, "res_var.tsv")
out <- data.frame(chr=chr,pos=p_names,sd=res_var)
write.table(out,f_out,col.names=F,row.names=F,sep="\t",quote=F)

########

meth_top <- get_top_meth(BSobj, v)

                                        # load genetic ancestry
f_ances <- here("inputs/genetic-ancestry/structure.out_ancestry_proportion_raceDemo_compare")
ances   <- get_ances(f_ances, meth_top$id_top)

                                        # residual after regressing out
                                        # Afr proportion
res <- res_ances(meth_top$meth, ances$gen_ances, ances$idx)

                                        # perform and plot pca
pc <- pca(res, output_path)
#PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10 
#0.05648 0.08545 0.11074 0.13106 0.15105 0.16687 0.17921 0.19071 0.20140 0.21163 

                                        # correlation of pc with ancestry
res <- corr_pc_ances(pc, ances$gen_ances, ances$idx, output_path)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
