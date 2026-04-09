## =============================================================================
## Load and Combine DLPFC WGBS Data
## =============================================================================
##
## PURPOSE:
##   This script loads per-chromosome BSseq objects for the dorsolateral
##   prefrontal cortex (DLPFC), combines them into a single genome-wide BSseq
##   object, and saves as an HDF5-backed SummarizedExperiment.
##
## INPUT:
##   - Per-chromosome HDF5SE from: inputs/wgbs-data/dlpfc/
##   - Phenotype data from: inputs/phenotypes/_m/phenotypes-AA.tsv
##
## OUTPUT:
##   - Combined BSseq object: inputs/wgbs-data/combined/_m/dlpfc_bsseq_h5se/
##   - Filtered phenotype data: inputs/wgbs-data/combined/_m/dlpfc_phenotypes.csv
##
## NOTE:
##   DLPFC data is already stored as HDF5SummarizedExperiment (not RDA),
##   so the loading method differs from caudate.
## =============================================================================
suppressMessages({
    library(here)
    library(bsseq)
    library(HDF5Array)
    library(DelayedArray)
    library(BiocParallel)
    library(SummarizedExperiment)
})

## Configuration
                                        # Set HDF5 options for efficient storage
options(HDF5Array.verbose = FALSE,
        HDF5Array.chunkdim = c(5000, 64),
        HDF5Array.max.block.size = 1e8)

region      <- "dlpfc"
output_path <- here("inputs/wgbs-data/combined/_m")

if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

chromosomes <- c(paste0("chr", 1:22), "chrX")
keep_assays <- c("M", "Cov")

## Load and Filter Phenotype Data
message("Loading phenotype data...")
pheno_file     <- here("inputs/phenotypes/_m", "phenotypes-all.tsv")
phenotype_data <- readr::read_tsv(pheno_file, show_col_types = FALSE) |>
    dplyr::filter(agedeath >= 17, region == !!region)

sample_ids <- phenotype_data$brnum
message(sprintf("Found %d samples for %s region", length(sample_ids), region))

## Define Helper Functions
harmonize_coldata_simple <- function(bs_list) {
    stopifnot(length(bs_list) > 0)

    common_cols <- Reduce(intersect, lapply(bs_list, \(b) colnames(colData(b))))
    if (length(common_cols) == 0) {
        stop("No common colData columns across objects.")
    }

    bs_list <- lapply(bs_list, \(b) {
        colData(b) <- colData(b)[, common_cols, drop = FALSE]; b
    })

    is_num_like <- function(x) inherits(x, c("numeric", "integer", "logical"))

    target_is_numeric <- vapply(common_cols, function(nm) {
        all(vapply(bs_list, \(b) is_num_like(colData(b)[[nm]]), logical(1)))
    }, logical(1))
    names(target_is_numeric) <- common_cols

    bs_list <- lapply(bs_list, function(b) {
        for (nm in common_cols) {
            if (target_is_numeric[[nm]]) {
                colData(b)[[nm]] <- as.numeric(colData(b)[[nm]])
            } else {
                colData(b)[[nm]] <- as.character(colData(b)[[nm]])
            }
        }
        b
    })
    return(bs_list)
}

## Load Per-Chromosome Data in Parallel
param <- MulticoreParam(
    workers = 6, progressbar = TRUE, stop.on.error = TRUE
)

message("Loading per-chromosome data in parallel...")
bs_list <- bplapply(chromosomes, function(chr) {
    message("  Loading ", chr, "...")
                                        # DLPFC uses HDF5SE format
    h5_dir <- here("inputs/wgbs-data", region, "_m",
                   paste0(region, "_", chr, "_BSobj"))
    BSobj  <- loadHDF5SummarizedExperiment(dir = h5_dir)

    if ("coef" %in% assayNames(BSobj)) { assays(BSobj)[["coef"]] <- NULL }
    assays(BSobj) <- assays(BSobj)[keep_assays]

                                        # Subset to samples of interest
    keep_cols <- colData(BSobj)$brnum %in% sample_ids
    BSobj     <- BSobj[, keep_cols, drop = FALSE]

                                        # Reorder samples to match sample_ids order
                                        # This ensures consistent sample order across chromosomes
    sample_order <- match(sample_ids, colData(BSobj)$brnum)
    sample_order <- sample_order[!is.na(sample_order)]
    BSobj        <- BSobj[, sample_order, drop = FALSE]
    BSobj
}, BPPARAM = param)
names(bs_list) <- chromosomes
gc()

## Harmonize and Combine Chromosomes
message("Harmonizing column data across chromosomes...")
bs_list <- harmonize_coldata_simple(bs_list)

message("Combining chromosomes into single object...")

grl <- do.call(c, lapply(bs_list, function(b) {
    rr <- rowRanges(b)
    mcols(rr) <- NULL
    rr
}))
grl <- as(grl, "GRangesList")

all_ranges <- unlist(grl)
all_assays <- lapply(keep_assays, function(nm) {
    do.call(DelayedArray::rbind, lapply(bs_list, \(b) assay(b, nm)))
})
names(all_assays) <- keep_assays

all_M   <- all_assays[["M"]]
all_Cov <- all_assays[["Cov"]]
rm(all_assays)

                                        # Get consistent sample names from first chromosome
                                        # (all chromosomes now have same sample order)
final_sample_names <- colData(bs_list[[1]])$brnum

                                        # Set column names to match expected sample names
colnames(all_M)   <- final_sample_names
colnames(all_Cov) <- final_sample_names

## Construct Final BSseq Object
message("Constructing combined BSseq object...")
bsseq_full <- BSseq(
    gr          = all_ranges,
    M           = all_M,
    Cov         = all_Cov,
    sampleNames = final_sample_names
)
colData(bsseq_full) <- colData(bs_list[[1]])
rownames(colData(bsseq_full)) <- colData(bsseq_full)$brnum

                                        # Clean up memory
rm(bs_list, all_ranges, all_M, all_Cov)
gc()

## Save Outputs
message("Saving filtered phenotype data...")
write.csv(
    phenotype_data, row.names = FALSE,
    file = file.path(output_path, paste0(region, "_phenotypes.csv"))
)

message("Saving BSseq as HDF5SummarizedExperiment...")
out_h5se <- file.path(output_path, paste0(region, "_bsseq_h5se"))
if (dir.exists(out_h5se)) unlink(out_h5se, recursive = TRUE)

saveHDF5SummarizedExperiment(
    bsseq_full, dir = out_h5se, replace = TRUE, verbose = FALSE
)

message(sprintf("Successfully saved %s BSseq object to: %s", region, out_h5se))

## Reproducibility Information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
