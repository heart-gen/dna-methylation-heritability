## =============================================================================
## Estimate cell type proportions using MuSiC on bulk RNA-seq data.
## =============================================================================
suppressPackageStartupMessages({
    library("here")
    library("MuSiC")
    library("dplyr")
    library("sessioninfo")
    library("SingleCellExperiment")
})

## Source shared utilities
source(here("inputs/cell_proportions/_h/deconvolution_utils.R"))
source(here("inputs/cell_proportions/_h/celltype_mapping.R"))

## Functions
run_music <- function(bulk_matrix, sce, marker_genes = NULL, verbose = FALSE) {
    if (verbose) cat("Run MuSiC deconvolution...\n")

                                        # Filter to common genes
    common_genes <- intersect(rownames(sce), colnames(bulk_matrix))
    if (!is.null(marker_genes)) {
        common_genes <- intersect(common_genes, marker_genes)
    }

    if (length(common_genes) < 50) {
        stop("Too few common genes between bulk and reference: ",
             length(common_genes))
    }

    if (verbose) cat("  Using", length(common_genes), "common genes\n")

                                        # Subset to common genes
    sce_sub  <- sce[common_genes, ]
    bulk_sub <- bulk_matrix[, common_genes]

                                        # Run MuSiC (MuSiC2 API: bulk.mtx + sc.sce)
    if (verbose) cat("  Running MuSiC estimation...\n")
    music_result <- tryCatch({
        music_prop(
            bulk.mtx = t(bulk_sub), sc.sce = sce_sub, clusters = "cell_type",
            samples = "donor", verbose = verbose
        )
    }, error = function(e) {
        message("MuSiC estimation failed: ", e$message)
        return(NULL)
    })

    if (is.null(music_result)) {
        stop("MuSiC deconvolution failed")
    }

                                        # Extract proportions
    result <- music_result$Est.prop.weighted

    if (verbose) cat("MuSiC deconvolution complete.\n")
    return(result)
}

load_bulk_data <- function(region, verbose = FALSE) {
    region_map <- list(
        "dlpfc" = "rse-gene.bsp2.dlpfc-n500.gencode-v47.RData",
        "hippocampus" = "rse-gene.bsp2.hippocampus-n452.gencode-v47.RData",
        "caudate" = "rse-gene.bsp3.caudate-n487.gencode-v47.RData"
    )

    if (!region %in% names(region_map)) {
        stop("Unknown region: ", region, ". Expected: caudate, dlpfc, or hippocampus")
    }

    if (verbose) cat("Loading bulk data for region:", region, "...\n")
    fn <- here("inputs/counts", region_map[[region]])
    load(fn, verbose = verbose)
    rse_df <- rse; rm(rse)
    return(rse_df)
}

bulk_rse_to_matrix <- function(rse_df, sce, verbose = FALSE) {
    if (verbose) cat("Converting bulk RSE to matrix...\n")
    bulk_counts <- assays(rse_df)$counts
    bulk_matrix <- t(as.matrix(bulk_counts))

    gene_names <- make.unique(rowData(rse_df)$gene_name)
    colnames(bulk_matrix) <- gene_names

    sample_names <- colnames(rse_df)
    if ("BrNum" %in% colnames(colData(rse_df))) {
        sample_names <- colData(rse_df)$BrNum
    }
    rownames(bulk_matrix) <- make.unique(as.character(sample_names))
    rownames(sce) <- rowData(sce)$gene_name
    return(bulk_matrix)
}

## Main Execution
args <- commandArgs(trailingOnly = TRUE)
region <- if (length(args) > 0) args[1] else "caudate"

                                        # Set seed for reproducibility
SEED <- 13
set.seed(SEED)

                                        # Load reference and bulk data
sce           <- load_sn_reference(region = region, verbose = TRUE)
bulk_rse      <- load_bulk_data(region, verbose = TRUE)
bulk_matrix   <- bulk_rse_to_matrix(bulk_rse, sce, verbose = TRUE)
rownames(sce) <- rowData(sce)$gene_name

                                        # Get marker genes
file_suffix  <- paste0("-", region)
marker_genes <- generate_marker_genes(sce, bulk_matrix, save_outputs = TRUE,
                                      file_suffix = region)
length(marker_genes)

                                        # Calculate cell proportions using MuSiC
est_prop_matrix <- run_music(bulk_matrix, sce, marker_genes = marker_genes,
                             verbose = TRUE)

                                        # Convert to data frame
music_prop <- est_prop_matrix |> as.data.frame() |>
    tibble::rownames_to_column("sample_id")

                                        # Save proportions
music_prop |>
    tidyr::pivot_longer(!sample_id, names_to = "cell_type",
                        values_to = "proportion") |>
    data.table::fwrite(paste0("music-proportions", file_suffix, ".tsv"),
                       sep = "\t")

## Session information
message("=== Reproducibility information ===")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
