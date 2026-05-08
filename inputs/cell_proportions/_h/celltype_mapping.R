## =============================================================================
## CELL TYPE MAPPING AND REFERENCE LOADING UTILITIES
## =============================================================================
suppressPackageStartupMessages({
    library("here")
    library("dplyr")
    library("SingleCellExperiment")
})

map_celltype_to_board <- function(cellType, region = "caudate", verbose = FALSE) {
    if (verbose) cat("Map cell types for region:", region, "\n")

                                        # Base mapping common to all regions
    mapped <- case_when(
        grepl("^drop", cellType) ~ NA_character_,
        grepl("^Astro", cellType) ~ "Astro",
        grepl("^Excit", cellType) ~ "Excit",
        grepl("^Inhib", cellType) ~ "Inhib",
        grepl("^Oligo", cellType) ~ "Oligo",
        cellType %in% c("Micro", "Micro_resting") ~ "Micro",
        cellType %in% c("Macrophage", "Tcell") ~ "Immune",
        cellType == "Mural" ~ "Mural",
        cellType %in% c("OPC", "OPC_COP") ~ "OPC",
        TRUE ~ NA_character_
    )

                                        # Region-specific mappings
    if (region == "caudate") {
        ## Nucleus Accumbens has MSN (medium spiny neurons)
        mapped <- case_when(
            !is.na(mapped) ~ mapped,
            grepl("^MSN.D1", cellType) ~ "D1-SPN",
            grepl("^MSN.D2", cellType) ~ "D2-SPN",
            TRUE ~ NA_character_
        )
    }

    return(mapped)
}

load_sn_reference <- function(region = "caudate", verbose = FALSE) {
    if (verbose) cat("Preparing data for region:", region, "\n")

                                        # Load region-specific SCE data
    sce <- switch(region,
        "caudate" = {
            load(here("inputs/sn-references/_m/SCE_NAc-n8_tran-etal.rda"))
            sce.nac.tran
        },
        "dlpfc" = {
            load(here("inputs/sn-references/_m/SCE_DLPFC-n3_tran-etal.rda"))
            sce.dlpfc.tran
        },
        "hippocampus" = {
            load(here("inputs/sn-references/_m/SCE_HPC-n3_tran-etal.rda"))
            sce.hpc.tran
        },
        stop("Unknown region: ", region, ". Expected: caudate, dlpfc, or hippocampus")
    )

                                        # Apply cell type mapping
    sce$cell_type <- map_celltype_to_board(sce$cellType, region = region)

                                        # Filter out excluded cells (NA cell types)
    sce <- sce[, !is.na(sce$cell_type)]

                                        # Convert to factor and drop unused levels
    sce$cell_type <- factor(sce$cell_type) |> droplevels()

    if (verbose) {
        cat("Cell type distribution:\n")
        print(table(sce$cell_type))
        cat("\nDonors:", length(unique(sce$donor)), "\n")
    }
    return(sce)
}
