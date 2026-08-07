#!/usr/bin/env Rscript
# LIBD-style covariates for AA-only tensorQTL (Level 3).
# mod = ~ Sex + Dx + Age + snpPC1 + snpPC2 + snpPC3
# + k expression PCs, k = sva::num.sv(TMM log2-CPM phenotype matrix, mod)

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(dplyr)
    library(sva)
    library(sessioninfo)
})

option_list <- list(
    make_option("--prepared-dir", dest = "prepared_dir", type = "character"),
    make_option("--outdir", type = "character", default = ""),
    make_option("--max-pcs", dest = "max_pcs", type = "integer", default = 50)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!nzchar(opt$prepared_dir)) stop("--prepared-dir is required")
prepared_dir <- opt$prepared_dir
outdir <- if (nzchar(opt$outdir)) opt$outdir else file.path(dirname(prepared_dir), "standard")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pheno <- data.table::fread(file.path(prepared_dir, "phenotypes.tsv")) |>
    dplyr::mutate(
        sample_id = as.character(BrNum),
        Sex = factor(as.character(Sex), levels = c("F", "M")),
        Dx = factor(as.character(Dx), levels = c("Control", "SCZD")),
        Age = as.numeric(Age),
        snpPC1 = as.numeric(snpPC1),
        snpPC2 = as.numeric(snpPC2),
        snpPC3 = as.numeric(snpPC3)
    )
if (any(is.na(pheno$Sex)) || any(is.na(pheno$Dx))) {
    stop("Sex/Dx has unexpected levels after factoring")
}

expr_dt <- data.table::fread(file.path(prepared_dir, "normalized_expression.tsv.gz"))
feat <- expr_dt[["feature_id"]]
expr <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr) <- feat
expr <- expr[, pheno$sample_id, drop = FALSE]

# Phenotype is already TMM log2-CPM
mod <- model.matrix(~ Sex + Dx + Age + snpPC1 + snpPC2 + snpPC3, data = pheno)

n_pc <- 0L
pcs <- NULL
pc_error <- NA_character_
tryCatch({
    if (nrow(expr) > 50000) {
        n_pc <- sva::num.sv(expr, mod, method = "be", vfilter = 50000)
    } else {
        n_pc <- sva::num.sv(expr, mod, method = "be")
    }
    n_pc <- min(as.integer(n_pc), as.integer(opt$max_pcs))
    if (n_pc > 0) {
        pca_df <- prcomp(t(expr), center = TRUE, scale. = FALSE)
        pcs <- as.data.frame(pca_df$x[, seq_len(n_pc), drop = FALSE])
        colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
    }
}, error = function(e) {
    pc_error <<- conditionMessage(e)
})

# Drop intercept and SexF reference column pattern: keep SexM, DxSCZD, Age, snpPCs
cov_fixed <- as.data.frame(mod[, colnames(mod) != "(Intercept)", drop = FALSE],
                           check.names = FALSE)
if (!is.null(pcs)) cov_fixed <- cbind(cov_fixed, pcs)

# TensorQTL-friendly: samples as rows
cov_out <- cbind(sample_id = pheno$sample_id, cov_fixed)

# Also write FastQTL-style transposed covariates
mat <- as.matrix(cov_out[, -1, drop = FALSE])
rownames(mat) <- cov_out$sample_id
cov_t <- data.table::as.data.table(t(mat), keep.rownames = "ID")

data.table::fwrite(cov_out, file.path(outdir, "covariates.txt"), sep = "\t", na = "NA")
data.table::fwrite(cov_t, file.path(outdir, "genes.combined_covariates.txt"), sep = "\t", na = "NA")

diag <- data.frame(
    region = NA_character_,
    feature = "genes",
    model = "standard_libd_aa",
    n_samples = nrow(pheno),
    n_features = nrow(expr),
    n_expr_pcs = n_pc,
    pc_error = pc_error,
    formula = "~ Sex + Dx + Age + snpPC1 + snpPC2 + snpPC3 + PC1..k",
    cohort_policy = "AA_AgeGt13_ControlSCZD_genotyped"
)
data.table::fwrite(diag, file.path(outdir, "covariate_diagnostics.tsv"), sep = "\t")
print(diag)
sessioninfo::session_info()
