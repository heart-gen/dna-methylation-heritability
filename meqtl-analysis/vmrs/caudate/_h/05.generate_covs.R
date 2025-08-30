## Prepare covariates for FastQTL / tensorQTL analysis
## Requires the 'sva' package
suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(SummarizedExperiment)
})

#### Functions
getRPKM <- function(rse, length_var = "bp_length", mapped_var = NULL) {
    mapped <- if (!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)
    bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)
    len <- if (!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))
    wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    return(assays(rse)$counts / (wid / 1000) / (bg / 1e6))
}

load_counts <- function(){
    caud8_file <- here("inputs/counts/_m",
                       "rse-gene.bsp3.caudate-n487.gencode-v47.RData")
    load(caud8_file)
    return(rse)
}

filter_data <- function(){
                                        # Load counts
    rse_df  <- load_counts()
                                        # Sample selection
    keepInd <- which((colData(rse_df)$Age > 13) &
                     (colData(rse_df)$dropped == "f") &
                     (colData(rse_df)$Race %in% c("AA", "CAUC")))
    rse_df  <- rse_df[,keepInd]
    rpkm_df <- getRPKM(rse_df, "length")
    rse_df  <- rse_df[rowMeans(rpkm_df) > 0.2,]; rm(rpkm_df)
    return(rse_df)
}
memDF <- memoise::memoise(filter_data)

normalize_data <- function(){
    rse_df <- memDF()
    norm_df <- getRPKM(rse_df, 'length')
    colnames(norm_df) <- gsub("\\_.*", "", colnames(norm_df))
    return(norm_df)
}

load_pca <- function(){
    gfile <- here("inputs/genotypes/_m/TOPMed_LIBD.eigenvec")
    pca   <- data.table::fread(gfile) |>
        rename_at(.vars = vars(starts_with("PC")),
                  function(x){sub("PC", "snpPC", x)}) |>
        mutate_if(is.character, as.factor)
    return(pca)
}

get_pheno <- function(){
    return(colData(memDF()) |> as.data.frame() |>
           inner_join(load_pca(), by=c("BrNum"="#FID"), multiple="all") |>
           distinct(RNum, .keep_all = TRUE) |> mutate(ids=RNum) |>
           tibble::column_to_rownames("ids"))
}

cal_pca <- function(norm_df, new_pd, mod){
    mat    <- log2(norm_df[, new_pd$RNum]+1)
    pca_df <- prcomp(t(mat))
    if(dim(norm_df)[1] > 50000){
        k  <- sva::num.sv(as.matrix(mat), mod, vfilter=50000)
    } else {
        k  <- sva::num.sv(as.matrix(mat), mod)
    }
    pcs    <- pca_df$x[,1:k]
    return(pcs)
}

#### Main

                                        # Load data
norm_df <- normalize_data()

                                        # Full model
new_pd  <- get_pheno()
mod     <- model.matrix(~Dx + Sex + Age + snpPC1 + snpPC2 + snpPC3,
                        data = new_pd)

                                        # PCA
pcs    <- cal_pca(norm_df, new_pd, mod)

                                        # Subset samples
sample_df <- data.table::fread("../_m/sample_id_to_brnum.tsv")

                                        # Extract covariates
covs <- cbind(mod[,c(-1, -2)], pcs) |> as.data.frame() |>
    tibble::rownames_to_column("RNum") |>
    inner_join(sample_df, by=c("RNum")) |> select(-"RNum") |>
    rename("BrNum"="ID") |> tibble::column_to_rownames("ID") |>
    t() |> as.data.frame() |> tibble::rownames_to_column("ID")
covs <- covs[, c("ID", sample_df$BrNum)]

                                        # Save file
data.table::fwrite(covs, "genes.combined_covariates.txt", sep='\t')

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
