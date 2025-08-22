## Script to prepare the data
suppressMessages({
    library(here)
    library(dplyr)
    library(SummarizedExperiment)
})

#### Phenotype data
                                        # Download from eQTL browser
caud8_file <- here("inputs/counts/_m",
                   "rse-gene.bsp3.caudate-n487.gencode-v47.RData")
load(caud8_file)
                                        # Select variables of interest
                                        # Update after confounding analysis
fields   <- c('BrNum', 'RNum', 'Region', 'RIN', 'Age', 'Sex', 'Race',
              'Dx', 'MoD', "dropped", "mito_mapping_rate", 'mapped_percent')
pheno    <- colData(rse)[,fields]
                                        # Basic filtering
pheno <- filter(as.data.frame(pheno), Age > 13, dropped == "f",
                Race %in% c("AA", "CAUC"))
                                        # Save file
as.data.frame(pheno) |>
    write.csv(file = "phenotypes.csv")

#### Normalized counts
keepIndex <- which(rse$Age > 13 & (rse$dropped == "f") &
                   rse$Race %in% c("AA", "CAUC"))
rse  <- rse[, keepIndex]
rm(keepIndex)
                                        # Generate DGE list
x      <- edgeR::DGEList(counts=assays(rse)$counts[, pheno$RNum],
                         genes=rowData(rse), samples=pheno)
                                        # Filter by expression
design <- model.matrix(~Dx, data=x$samples)
keep.x <- edgeR::filterByExpr(x, design=design)
print(paste('There are:', sum(keep.x), 'features left!', sep=' '))
x      <- x[keep.x, , keep.lib.sizes=FALSE]
                                        # Normalize library size
x      <- edgeR::calcNormFactors(x, method="TMM")
                                        # Normalize counts
cpm    <- edgeR::cpm(x)
                                        # Save data
as.data.frame(cpm) |> tibble::rownames_to_column("feature_id") |>
    write.table(file=gzfile("normalized_expression.txt.gz"),
                sep="\t", quote=FALSE, row.names=FALSE)
rm(keep.x, design, cpm)

#### Export annotation
genes_to_keep <- rownames(x$genes)

annot <- rowRanges(rse) |> as.data.frame() |> distinct() |>
    tibble::rownames_to_column("feature_id") |>
    filter(gene_id %in% genes_to_keep) |>
    select(seqnames, start, end, feature_id, strand) |>
    mutate(index=feature_id) |> tibble::column_to_rownames("index")
                                        # Save data
data.table::fwrite(annot, "feature.bed", sep='\t', row.names=TRUE)

#### R session information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
