##### Performs GO enrichment analysis for all brain regions #####
suppressPackageStartupMessages({
    library(rGREAT)
    library(here)
    library(dplyr)
    library(plyranges)
    library(purrr)
    library(KEGGREST)
    library(reactome.db)
    library(UniProtKeywords)
})

## --- Function --- ##
filter_heritability <- function(tissue,
                                heritability_filter = "all",
                                pval_filter = NULL,
                                apply_h2_filter = FALSE) {
                                        # Read in summary table
    vmr_test <- "/projects/p32505/users/alexis/projects/dna-methylation-heritability/heritability/gcta/caudate/_m/summary_to_del/greml_summary.tsv"
    vmr_file <- here("heritability/gcta", tissue,
                     "_m/summary/greml_summary.tsv")
    vmr <- read.table(vmr_test, sep = "\t", header = TRUE)

                                        # Apply FDR correction
    vmr$Pval_Variance <- p.adjust(vmr$Pval_Variance, method = "fdr")

                                        # Filter by heritability
    if (apply_h2_filter) {
        h2_thresh  <- 0.75
        if (heritability_filter == "heritable") {
            vmr <- dplyr::filter(vmr, Sum.of.V.G._Vp_Variance >= h2_thresh)
        } else if (heritability_filter == "non_heritable") {
            vmr <- dplyr::filter(vmr, Sum.of.V.G._Vp_Variance < h2_thresh)
        }
    }

                                        # Filter by p-value
    if (pval_filter != "NULL") {
        if (heritability_filter == "heritable") {
            vmr <- dplyr::filter(vmr, Pval_Variance < pval_filter)
        } else if (heritability_filter == "non_heritable") {
            vmr <- dplyr::filter(vmr, Pval_Variance > pval_filter)
        }
    }

    return(vmr)
}

## --- GO enrichment --- ##
load_vmr_background <- function(tissue) {
                                        # Load the regions tested as background
    vmr_file <- here("heritability/gcta", tissue, "_m/vmr_list.txt")
    vmr_df   <- read.table(vmr_file)
    colnames(vmr_df) <- c("seqnames", "start", "end")
    vmr_gr <- plyranges::as_granges(vmr_df)
    seqlevels(vmr_gr) <- paste0("chr", seqlevels(vmr_gr))
    return(vmr_gr)
}

get_enrichment <- function(vmr_filtered, tissue, filter_label) {
    colnames(vmr_filtered) <- c("seqnames", "start", "end")
    vmr <- plyranges::as_granges(vmr_filtered)
    seqlevels(vmr) <- paste0("chr", seqlevels(vmr))

    gene_sets <- list(
        "GO:BP" = "GO:BP",
        "GO:MF" = "GO:MF",
        KEGG = split(gsub("hsa:", "", names(keggLink("pathway","hsa"))),
                     gsub("path:", "", keggLink("pathway", "hsa"))),
        reactome = as.list(reactomePATHID2EXTID),
        uniprot = load_keyword_genesets(9606),
        "msigdb:C7:IMMUNESIGDB" = "msigdb:C7:IMMUNESIGDB"
    )

    background_df <- load_vmr_background(tissue)
    for (gs in names(gene_sets)) {
      message("Running GREAT for ", gs)
      res <- great(vmr, gene_sets[[gs]], "RefSeq:hg38",
                   background = background_df)
      tb  <- getEnrichmentTable(res)
      new_gs <- gsub(":", "_", gs)
      outfile <- paste0(tissue, "_", filter_label, "_", new_gs, ".csv")
      write.csv(
        tb,
        file = here("heritability/gcta/tissue_comparison",
                    "functional_enrichment/_m", outfile),
        row.names = FALSE
      )
    }
}

# Main
tissues              <- c("caudate", "dlpfc", "hippocampus")
heritability_filters <- c("all", "heritable", "non_heritable")
pval_filters         <- c("NULL", 0.05, 0.1, 0.25)
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
            parts <- c(parts, "pval", as.character(pval))
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
      message("No data left after filtering for ", tissue, " - ",label)
    }
  }
})

# Reproducibility info
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
