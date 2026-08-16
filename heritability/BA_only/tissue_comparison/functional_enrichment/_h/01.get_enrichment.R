##### Performs GO enrichment analysis for all brain regions #####
suppressPackageStartupMessages({
    library(rGREAT)
    library(dplyr)
    library(plyranges)
    library(here)
    library(purrr)
    library(KEGGREST)
    library(reactome.db)
    library(UniProtKeywords)
})

## --- Function --- ##
filter_heritability <- function(tissue, heritability_filter) {
  
                                        # Read in summary table
    vmr_file <- here("heritability/elastic_net_model/BA_only",
                      tissue, "_m", paste0(tissue, "_summary_elastic-net.tsv"))
    vmr <- read.table(vmr_file, sep = "\t", header = TRUE) %>% na.omit()

                                        # Filter by h2_category
    h2_thresh  <- 0.1
    r2_thresh  <- 0.3
    if (heritability_filter == "heritable") {
        vmr <- vmr %>% filter(h2_unscaled >= h2_thresh & r_squared_cv > r2_thresh)
    } else if (heritability_filter == "non_heritable") {
        vmr <- vmr %>% filter(h2_unscaled < h2_thresh & r_squared_cv > r2_thresh)
    } else if (heritability_filter == "low_prediction") {
        vmr <- vmr %>% filter(r_squared_cv <= r2_thresh)
    } else if (heritability_filter == "all") {
        vmr <- vmr
    }
  
    return(vmr)
}

## --- GO enrichment --- ##
load_vmr_background <- function(tissue) {
                                        # Load the regions tested as background
    vmr_file <- here("vmr-analysis", tissue, "_m/vmr.bed")
    vmr_df   <- read.table(vmr_file)
    colnames(vmr_df) <- c("seqnames", "start", "end")
    vmr_gr <- plyranges::as_granges(vmr_df)
    seqlevels(vmr_gr) <- paste0("chr", seqlevels(vmr_gr))
    return(vmr_gr)
}

get_enrichment <- function(vmr_filtered, tissue, hfilter) {
    vmr_filtered <- vmr_filtered[, c("chrom", "start", "end")]
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
      out_path <- here("heritability", "elastic_net_model", "BA_only", 
                      "tissue_comparison", "functional_enrichment", 
                      "_m", new_gs)
      if (!dir.exists(out_path)) {
        dir.create(out_path, recursive = TRUE)
      }
      outfile <- paste0(tissue, "_", hfilter, ".csv")
      write.csv(tb, file = file.path(out_path, outfile), row.names = FALSE)
    }
}

# Main
tissues              <- c("caudate", "dlpfc", "hippocampus")
heritability_filters <- c("all", "heritable", "non_heritable", "low_prediction")

# Run analysis
for (tissue in tissues) {
  for (hfilter in heritability_filters) {
    message("Running enrichment: ", tissue, " - ", hfilter)
      
    # Stratify VMRs based on heritability 
    vmr_filtered <- filter_heritability(tissue, hfilter)
    # Get enrichment for remaining data after filtering
    if (nrow(vmr_filtered) > 0) {
      get_enrichment(vmr_filtered, tissue, hfilter)
    } else {
      message("No data left after filtering for ", tissue, " - ", hfilter)
    }
  }
}

# Reproducibility info
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
