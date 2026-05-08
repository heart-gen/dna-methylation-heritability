## =============================================================================
## SHARED DECONVOLUTION UTILITIES
## =============================================================================
suppressPackageStartupMessages({
    library("dplyr")
    library("DeconvoBuddies")
    library("SingleCellExperiment")
})

generate_marker_genes <- function(sce, bulk_matrix, verbose = FALSE,
                                  save_outputs = FALSE, file_suffix = "dlpfc") {
    if (verbose) cat("Select marker genes...\n")
    rownames(sce) <- rowData(sce)$gene_name
                                        # Creating mean_ratios
    ratios <- get_mean_ratio(sce, cellType_col = "cell_type",
                             assay_name = "logcounts",
                             gene_ensembl = "gene_id",
                             gene_name = "gene_name")
    if (save_outputs) {
        csv_fn <- paste0("marker_stats_genes.", file_suffix, ".csv")
        pdf_fn <- paste0("top2-marker-genes.", file_suffix, ".pdf")

        write.csv(ratios, file = csv_fn, row.names = FALSE)

        plot_marker_express_ALL(
            sce, ratios, n_genes = 2, pdf_fn = pdf_fn,
            cellType_col = "cell_type",
            color_pal = NULL, plot_points = FALSE
        )
    }

    ## creating marker_list of top 25 genes
    marker_genes <- ratios |>
        filter(MeanRatio.rank <= 25, gene_name %in% colnames(bulk_matrix)) |>
        pull(gene_name)
    return(marker_genes)
}
