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

pca <- function(meth, output_path) {
    pc <- prcomp(t(meth), scale. = TRUE, center = TRUE)
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

## Main
                                        # create output directories if they  
                                        # don't exist
output_path <- here("covar-analysis", "caudate", "no_snp_pc", "_m", "pca", paste0("chr_", chr))
if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

                                        # load sd of raw DNAm
load(here("covar-analysis", "caudate", "no_snp_pc", "_m", "cpg", paste0("chr_", chr), "stats.rda"))
v <- data.frame(chr = chr, start = start(BSobj), sd = sds)

                                        # get top 1M variable CpG
meth_top <- get_top_meth(BSobj, v)

                                        # perform and plot pca
pc <- pca(meth_top$meth, output_path)
#    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10 
#0.06703 0.09666 0.11962 0.14005 0.15736 0.17040 0.18309 0.19547 0.20673 0.21757

                                        # correlation of pc with ancestry
#corr_pc_ances(pc, ances$gen_ances, ances$idx, output_path)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
