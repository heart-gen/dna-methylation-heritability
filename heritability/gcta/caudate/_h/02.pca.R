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
