#!/usr/bin/env Rscript
# Prepare AA-only BrainSeq gene expression for LIBD-style tensorQTL (Level 3).
# Sample filter: Age > 13, Race == AA, Dx %in% Control/SCZD, genotyped.
# Feature filter: edgeR::filterByExpr; phenotype: TMM log2-CPM.

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(dplyr)
    library(SummarizedExperiment)
    library(edgeR)
    library(sessioninfo)
})

option_list <- list(
    make_option("--region", type = "character", default = "caudate"),
    make_option("--outdir", type = "character"),
    make_option("--rse", type = "character", default = ""),
    make_option("--annotation", type = "character", default = ""),
    make_option("--genotype-psam", type = "character", default = ""),
    make_option("--eigenvec", type = "character", default = ""),
    make_option("--min-age", type = "double", default = 13),
    make_option("--prior-count", dest = "prior_count", type = "double", default = 0.25)
)
opt <- parse_args(OptionParser(option_list = option_list))
region <- tolower(opt$region)
if (!nzchar(opt$outdir)) stop("--outdir is required")
outdir <- opt$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

root <- "/projects/b1213/users/kynon/projects/dna-methylation-heritability"
rse_default <- list(
    caudate = file.path(root, "inputs/counts/rse-gene.bsp3.caudate-n487.gencode-v47.RData"),
    dlpfc = file.path(root, "inputs/counts/rse-gene.bsp2.dlpfc-n500.gencode-v47.RData"),
    hippocampus = file.path(root, "inputs/counts/rse-gene.bsp2.hippocampus-n452.gencode-v47.RData")
)
annot_default <- file.path(root, "inputs/counts/gene-annotation.tsv")
psam_default <- file.path(root, "inputs/genotypes/TOPMed_LIBD.AA.psam")
eigen_default <- file.path(root, "inputs/genotypes/TOPMed_LIBD.AA.eigenvec")

rse_path <- if (nzchar(opt$rse)) opt$rse else rse_default[[region]]
annot_path <- if (nzchar(opt$annotation)) opt$annotation else annot_default
psam_path <- if (nzchar(opt$`genotype-psam`)) opt$`genotype-psam` else psam_default
eigen_path <- if (nzchar(opt$eigenvec)) opt$eigenvec else eigen_default
if (is.null(rse_path) || !file.exists(rse_path)) stop("Missing RSE: ", rse_path)

load(rse_path)
if (!exists("rse")) stop("Expected object 'rse' in ", rse_path)
cd <- as.data.frame(colData(rse))
cd$RNum <- as.character(cd$RNum)
cd$BrNum <- as.character(cd$BrNum)
cd$Age <- as.numeric(cd$Age)
cd$Race <- as.character(cd$Race)
cd$Dx <- as.character(cd$Dx)
cd$Sex <- as.character(cd$Sex)

# Genotype IDs (psam may be headerless: BrNum, chip, SEX)
psam_head <- readLines(psam_path, n = 1)
if (grepl("FID|IID|BrNum", psam_head, ignore.case = TRUE)) {
    psam_raw <- data.table::fread(psam_path, header = TRUE)
} else {
    psam_raw <- data.table::fread(psam_path, header = FALSE)
}
geno_ids <- unique(as.character(psam_raw[[1]]))
geno_ids <- geno_ids[grepl("^Br", geno_ids)]

keep <- which(
    cd$Age > opt$`min-age` &
        cd$Race == "AA" &
        cd$Dx %in% c("Control", "SCZD") &
        !is.na(cd$BrNum) &
        cd$BrNum %in% geno_ids &
        !is.na(cd$Sex) &
        cd$RNum %in% colnames(rse)
)
if (length(keep) < 20) stop("Too few samples after AA-only filter: ", length(keep))
rse <- rse[, keep]
cd <- cd[keep, , drop = FALSE]

cd$Sex <- dplyr::case_when(
    cd$Sex %in% c("M", "Male", "1") ~ "M",
    cd$Sex %in% c("F", "Female", "2") ~ "F",
    TRUE ~ cd$Sex
)
cd$Dx <- factor(cd$Dx, levels = c("Control", "SCZD"))
cd$Sex <- factor(cd$Sex, levels = c("F", "M"))

counts <- as.matrix(assays(rse)$counts)
storage.mode(counts) <- "double"
colnames(counts) <- cd$BrNum
# Collapse rare duplicate BrNums by keeping first
if (any(duplicated(colnames(counts)))) {
    keep_u <- !duplicated(colnames(counts))
    counts <- counts[, keep_u, drop = FALSE]
    cd <- cd[keep_u, , drop = FALSE]
}

x <- edgeR::DGEList(counts = counts, samples = cd)
design <- model.matrix(~ Dx + Sex + Age, data = cd)
keep_feat <- edgeR::filterByExpr(x, design = design)
n_input <- nrow(x)
x <- edgeR::calcNormFactors(x[keep_feat, , keep.lib.sizes = FALSE], method = "TMM")
expr <- edgeR::cpm(x, log = TRUE, prior.count = opt$prior_count)
storage.mode(expr) <- "double"

# Ancestry PCs: join BrNum to eigenvec #FID
eigen <- data.table::fread(eigen_path)
fid_col <- if ("#FID" %in% names(eigen)) "#FID" else names(eigen)[1]
pc_cols <- grep("^PC[0-9]+$", names(eigen), value = TRUE)
if (length(pc_cols) < 3) stop("Need at least PC1-PC3 in ", eigen_path)
eigen <- eigen |>
    dplyr::mutate(BrNum = as.character(.data[[fid_col]])) |>
    dplyr::select(BrNum, dplyr::all_of(pc_cols))
names(eigen) <- sub("^PC", "snpPC", names(eigen))

pheno <- cd |>
    dplyr::mutate(SAMPLE_ID = RNum, BrNum = as.character(BrNum)) |>
    dplyr::left_join(eigen, by = "BrNum")
if (any(is.na(pheno$snpPC1))) {
    stop("Missing snpPCs for ", sum(is.na(pheno$snpPC1)), " samples")
}
pheno <- pheno[match(colnames(expr), pheno$BrNum), , drop = FALSE]

# Annotation / feature BED (gene body; TSS collapse done in BED maker)
annot <- data.table::fread(annot_path)
id_col <- if ("gene_id" %in% names(annot)) "gene_id" else names(annot)[1]
annot <- annot |>
    dplyr::filter(.data[[id_col]] %in% rownames(expr)) |>
    dplyr::mutate(
        feature_id = .data[[id_col]],
        chrom = as.character(chrom),
        chrom = ifelse(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))
    ) |>
    dplyr::select(chrom, start, end, feature_id, strand, dplyr::everything())

common <- intersect(rownames(expr), annot$feature_id)
expr <- expr[common, , drop = FALSE]
annot <- annot[match(common, annot$feature_id), , drop = FALSE]
pheno <- pheno[match(colnames(expr), pheno$BrNum), , drop = FALSE]

expr_df <- data.table::as.data.table(expr, keep.rownames = "feature_id")
data.table::fwrite(expr_df, file.path(outdir, "normalized_expression.tsv.gz"), sep = "\t")
data.table::fwrite(annot, file.path(outdir, "feature.bed"), sep = "\t")
data.table::fwrite(pheno, file.path(outdir, "phenotypes.tsv"), sep = "\t", na = "NA")
data.table::fwrite(
    pheno[, c("RNum", "BrNum")],
    file.path(outdir, "sample_id_to_brnum.tsv"),
    sep = "\t"
)
data.table::fwrite(
    data.frame(chr = paste0("chr", c(1:22, "X"))),
    file.path(outdir, "vcf_chr_list.txt"),
    sep = "\t",
    col.names = FALSE
)

summary <- data.frame(
    region = region,
    feature = "genes",
    cohort_policy = "AA_AgeGt13_ControlSCZD_genotyped",
    n_samples = ncol(expr),
    n_features_input = n_input,
    n_features = nrow(expr),
    min_age = opt$`min-age`,
    feature_filter = "edgeR::filterByExpr(~ Dx + Sex + Age)",
    norm_method = "TMM",
    phenotype = "log2-CPM",
    prior_count = opt$prior_count,
    genotype_prefix = "inputs/genotypes/TOPMed_LIBD.AA",
    rse_path = rse_path
)
data.table::fwrite(summary, file.path(outdir, "prep_summary.tsv"), sep = "\t")
print(summary)
sessioninfo::session_info()
