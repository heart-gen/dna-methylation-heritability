## =============================================================================
## Load and Combine Caudate WGBS Data
## =============================================================================
##
## PURPOSE:
##   This script loads per-chromosome BSseq objects for the caudate nucleus,
##   combines them into a single genome-wide BSseq object, and saves as an
##   HDF5-backed SummarizedExperiment for efficient downstream analysis.
##
## INPUT:
##   - Per-chromosome BSseq RDA files from: inputs/wgbs-data/caudate/
##   - Phenotype data from: inputs/phenotypes/_m/phenotypes-AA.tsv
##
## OUTPUT:
##   - Combined BSseq object: inputs/wgbs-data/combined/_m/caudate_bsseq_h5se/
##   - Filtered phenotype data: inputs/wgbs-data/combined/_m/caudate_phenotypes.csv
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
                                        # These settings control how data is
                                        # chunked and processed
options(HDF5Array.verbose = FALSE,
        HDF5Array.chunkdim = c(5000, 64),
        HDF5Array.max.block.size = 1e8)

                                        # Define region and output paths
region      <- "caudate"
output_path <- here("inputs/wgbs-data/combined/_m")

                                        # Create output directory if it doesn't exist
if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

                                        # Define chromosomes to process
chromosomes <- c(paste0("chr", 1:22), "chrX")

                                        # Select only the assays we need for analysis:
                                        # - M: Methylated read counts
                                        # - Cov: Total coverage (reads covering each CpG)
keep_assays <- c("M", "Cov")

## Load and Filter Phenotype Data
message("Loading phenotype data...")
pheno_file     <- here("inputs/phenotypes/_m", "phenotypes-all.tsv")
phenotype_data <- readr::read_tsv(pheno_file, show_col_types = FALSE) |>
    dplyr::filter(agedeath >= 17, region == !!region)

                                        # Get sample IDs to keep
sample_ids <- phenotype_data$brnum
message(sprintf("Found %d samples for %s region", length(sample_ids), region))

## Define Helper Functions
load_bsobj_safely <- function(file) {
    env <- new.env(parent = emptyenv())
    on.exit(rm(env), add = TRUE)
    load(file, envir = env)

    if (!exists("BSobj", envir = env, inherits = FALSE)) {
        stop("Object `BSobj` not found in: ", file)
    }
    return(env$BSobj)
}

harmonize_coldata_simple <- function(bs_list) {
    stopifnot(length(bs_list) > 0)

                                        # Find columns that exist in ALL objects
    common_cols <- Reduce(intersect, lapply(bs_list, \(b) colnames(colData(b))))
    if (length(common_cols) == 0) {
        stop("No common colData columns across objects.")
    }

                                        # Restrict colData to common columns only
    bs_list <- lapply(bs_list, \(b) {
        colData(b) <- colData(b)[, common_cols, drop = FALSE]; b
    })

                                        # Helper to check if a column is numeric-like
    is_num_like <- function(x) inherits(x, c("numeric", "integer", "logical"))

                                        # Determine which columns should be numeric
    target_is_numeric <- vapply(common_cols, function(nm) {
        all(vapply(bs_list, \(b) is_num_like(colData(b)[[nm]]), logical(1)))
    }, logical(1))
    names(target_is_numeric) <- common_cols

                                        # Coerce columns to consistent types
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
                                        # Construct path to the chromosome-specific RDA file
    rda_file <- here("inputs/wgbs-data", region, "_m",
                 paste0(region, "_", chr, "_BSobj.rda"))

    BSobj <- load_bsobj_safely(rda_file)

                                        # Remove the 'coef' assay if present (not needed)
    an <- assayNames(BSobj)
    if ("coef" %in% an) { assays(BSobj)[["coef"]] <- NULL }
    assays(BSobj) <- assays(BSobj)[keep_assays]

                                        # Subset to only our samples of interest
    keep_cols <- colData(BSobj)$brnum %in% sample_ids
    BSobj     <- BSobj[, keep_cols, drop = FALSE]
    BSobj
}, BPPARAM = param)
names(bs_list) <- chromosomes
gc()

## Harmonize and Combine Chromosomes
message("Harmonizing column data across chromosomes...")
bs_list <- harmonize_coldata_simple(bs_list)

message("Combining chromosomes into single object...")

                                        # Combine all genomic ranges
grl <- do.call(c, lapply(bs_list, function(b) {
    rr <- rowRanges(b)
    mcols(rr) <- NULL
    rr
}))
grl <- as(grl, "GRangesList")

all_ranges <- unlist(grl)

                                        # Combine assays by row-binding
all_assays <- lapply(keep_assays, function(nm) {
    do.call(DelayedArray::rbind, lapply(bs_list, \(b) assay(b, nm)))
})
names(all_assays) <- keep_assays

all_M   <- all_assays[["M"]]
all_Cov <- all_assays[["Cov"]]
rm(all_assays)

## Construct Final BSseq Object
message("Constructing combined BSseq object...")
bsseq_full <- BSseq(
    gr          = all_ranges,
    M           = all_M,
    Cov         = all_Cov,
    sampleNames = colData(bs_list[[1]])$brnum
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
