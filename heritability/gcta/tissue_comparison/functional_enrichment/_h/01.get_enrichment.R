##### Performs GO enrichment analysis for all brain regions #####
suppressPackageStartupMessages({
    library(rGREAT)
    library(here)
    library(dplyr)
    library(purrr)
})

# Function
filter_heritability <- function(tissue,
                                heritability_filter = "all",
                                pval_filter = NULL,
                                apply_h2_filter = FALSE) {
    # Read in summary table
    vmr <- read.table(here("heritability/gcta/",
                           paste0(tissue, "/_m/summary/greml_summary.tsv")),
                     sep = "\t", header = TRUE)

    # Filter by heritability
    if (apply_h2_filter) {
        top_10  <- quantile(vmr$Sum.of.V.G._Vp_Variance, 0.90, na.rm = TRUE)
        if (heritability_filter == "heritable") {
            vmr <- vmr %>% dplyr::filter(Sum.of.V.G._Vp_Variance >= top_10)
        } else if (heritability_filter == "non_heritable") {
            vmr <- vmr %>% dplyr::filter(Sum.of.V.G._Vp_Variance < top_10)
        }
    }

    # Filter by p-value
    if (pval_filter != "NULL") {
        if (heritability_filter == "heritable") {
            vmr <- vmr %>% dplyr::filter(Pval_Variance < pval_filter)
        } else if (heritability_filter == "non_heritable") {
            vmr <- vmr %>% dplyr::filter(Pval_Variance > pval_filter)
        }
    }

    return(vmr)
}

# GO enrichment
get_enrichment <- function(vmr_filtered, tissue, filter_label) {
    colnames(vmr_filtered) <- c("seqnames", "start", "end")
    vmr <- plyranges::as_granges(vmr_filtered)
    res <- great(vmr, "GO:BP", "RefSeq:hg38", background = as.character(1:22))
    tb  <- getEnrichmentTable(res)
    write.csv(
        tb,
        file = here("heritability/gcta/tissue_comparison/functional_enrichment/_m",
                    paste0(tissue, "_", filter_label, "_GO.csv")),
        row.names = FALSE
    )
}

# Main
tissues              <- c("caudate", "dlpfc", "hippocampus")
heritability_filters <- c("all", "heritable", "non_heritable")
pval_filters         <- c("NULL", 0.05)
apply_h2_options     <- c(TRUE, FALSE)

# Create all possible combinations
filter_grid <- expand.grid(
    hfilter          = heritability_filters,
    pval             = pval_filters,
    h2               = apply_h2_options,
    stringsAsFactors = FALSE
) %>%
    # Filter out invalid combinations
    dplyr::filter(
        !(hfilter == "all" & (pval != "NULL" | h2 != FALSE)),
        !(hfilter %in% c("heritable", "non_heritable") & pval == "NULL" & h2 == FALSE)
    ) %>%
    # Get labels for each combination
    rowwise() %>%
    mutate(
        label = {
            parts <- c(hfilter)
        if (h2) {
            parts <- c(parts, "h2")
        }
        if (pval != "NULL") {
            parts <- c(parts, "pval")
        }
        paste(parts, collapse = "_")
      }
    ) %>%
    ungroup()

# Run analysis
purrr::pwalk(filter_grid, function(hfilter, pval, h2, label) {
  for (tissue in tissues) {
    message("Running enrichment: ", tissue, " - ", label)
  
    # Apply the heritability and p-value filter
    vmr_filtered <- filter_heritability(tissue, hfilter, pval, h2)
    # Get enrichment for remaining data after filtering
    if (nrow(vmr_filtered) > 0) {
      get_enrichment(vmr_filtered, tissue, label)
    } else {
      message("No data left after filtering for ", tissue, " - ", label)
    }
  }
})

# Reproducibility info
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()