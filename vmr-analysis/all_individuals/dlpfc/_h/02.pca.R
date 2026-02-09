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

get_snp_pcs <- function(pc_file_path) {
    pc_cols <- paste0("snpPC", 1:10)
    snp_pcs <- fread(pc_file_path, header=TRUE)[
      agedeath >= 17 & region == "dlpfc",
      c("brnum", pc_cols), with = FALSE
    ]
    na_samples <- snp_pcs$brnum[rowSums(is.na(snp_pcs)) > 0]

    if(length(na_samples) > 0) {
      message("Removing samples with NA values: ", paste(na_samples))
      snp_pcs <- snp_pcs[!brnum %in% na_samples]
    }
    return(snp_pcs)
}

res_snp_pcs <- function(meth, snp_pcs) {
    meth <- meth[, colnames(meth) %in% snp_pcs$brnum]
    snp_pcs <- snp_pcs[match(colnames(meth), snp_pcs$brnum), 
                       paste0("snpPC", 1:3)]
    design <- cbind(1, as.matrix(snp_pcs))
    fit    <- limma::lmFit(meth, design)
    resids <- limma::residuals.MArrayLM(fit, meth)

    return(resids)
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

## Main
                                        # create output directories if they  
                                        # don't exist
output_path <- here("vmr-analysis", "all_individuals", "dlpfc",
                    "_m", "pca", paste0("chr_", chr))
if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

                                        # load sd of raw DNAm
load(here("vmr-analysis", "all_individuals", "dlpfc",
          "_m", "cpg", paste0("chr_", chr), "stats.rda"))
v <- data.frame(chr = chr, start = start(BSobj), sd = sds)

                                        # get top 1M variable CpG
meth_top <- get_top_meth(BSobj, v)

                                        # get SNP PCs
pc_file_path <- here("inputs/phenotypes/_m/phenotypes-all.tsv")
snp_pcs   <- get_snp_pcs(pc_file_path)

                                        # residual after regressing out
                                        # top 3 SNP PCs
res <- res_snp_pcs(meth_top$meth, snp_pcs)

                                        # perform and plot pca
pc <- pca(res, output_path)
#PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10 
#0.05648 0.08545 0.11074 0.13106 0.15105 0.16687 0.17921 0.19071 0.20140 0.21163 

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
