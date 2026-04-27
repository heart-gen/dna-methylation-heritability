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
filter_sites <- function(enet) {
  vmr <- na.omit(enet)
  vmr <- vmr %>%
    mutate(h2_category = case_when(
      r_squared_cv <= 0.3 ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.3 ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.3 ~ "Heritable"
    ),
    h2_category = factor(h2_category, levels = c("Heritable", 
                                                 "Non-heritable", 
                                                 "Low prediction"))
    )
  return(vmr)
}

filter_heritability <- function(vmr, heritability_filter) {
                                        # Filter by h2_category
    if (heritability_filter == "heritable") {
        vmr <- vmr %>% filter(h2_category == "Heritable")
    } else if (heritability_filter == "non_heritable") {
        vmr <- vmr %>% filter(h2_category == "Non-heritable")
    } else if (heritability_filter == "low_prediction") {
        vmr <- vmr %>% filter(h2_category == "Low prediction")
    } else if (heritability_filter == "all") {
        vmr <- vmr
    }
  
    return(vmr)
}

## --- GO enrichment --- ##
load_vmr_background <- function(enet, tissue, pop) {
                                        # Load the regions tested as background
    vmr_file <- here("vmr-analysis/all_individuals", tissue, "_m/vmr.bed")
    vmr_df   <- read.table(vmr_file)
    colnames(vmr_df) <- c("seqnames", "start", "end")

    if (pop == "matched") {
      enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                      paste0(tissue, "/_m/", tissue, "_summary_elastic-net_matched_r2_0.3.tsv"))
      shared_vmrs <- read.table(enet_file, sep = "\t", header = TRUE) %>% na.omit()
                                          # Keep only shared vmrs as background
      vmr_df <- vmr_df %>%
          semi_join(shared_vmrs, by = c("seqnames" = "chrom", "start", "end"))
    }

    vmr_gr <- plyranges::as_granges(vmr_df)
    seqlevels(vmr_gr) <- paste0("chr", seqlevels(vmr_gr))
    return(vmr_gr)
}

get_enrichment <- function(vmr_filtered, tissue, pop, hfilter) {
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

    background_df <- load_vmr_background(tissue, pop)
    for (gs in names(gene_sets)) {
      message("Running GREAT for ", gs)
      res <- great(vmr, gene_sets[[gs]], "RefSeq:hg38",
                   background = background_df)
      tb  <- getEnrichmentTable(res)
      new_gs <- gsub(":", "_", gs)
      out_path <- here("heritability", "elastic_net_model", "all_individuals", 
                      "tissue_comparison", "functional_enrichment", 
                      "_m", new_gs)
      if (!dir.exists(out_path)) {
        dir.create(out_path, recursive = TRUE)
      }
      outfile <- paste0(tissue, "_", hfilter, "_", pop, ".csv")
      write.csv(tb, file = file.path(out_path, outfile), row.names = FALSE)
    }
}

# Main
tissues              <- c("caudate", "dlpfc", "hippocampus")
populations          <- c("AA", "EA", "matched")
heritability_filters <- c("all", "heritable", "non_heritable", "low_prediction")

# Run analysis
for (pop in populations) {
  for (tissue in tissues) {

    # Read in summary table
    if (pop %in% c("AA", "EA")) {
      enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                        paste0(tissue, "/_m/", tissue, "_summary_elastic-net_", pop, ".tsv"))
      enet <- read.table(enet_file, sep = "\t", header = TRUE)
      vmr  <- filter_sites(enet)

    } else if (pop == "matched") {
      enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                        paste0(tissue, "/_m/", tissue, "_summary_elastic-net_matched_r2_0.3.tsv"))
      vmr <- read.table(enet_file, sep = "\t", header = TRUE) %>% na.omit()
    }
    
    for (hfilter in heritability_filters) {
      message("Running enrichment: ", pop, " - ", tissue, " - ", hfilter)
        
      # Stratify VMRs based on heritability 
      vmr_filtered <- filter_heritability(vmr, hfilter)
      # Get enrichment for remaining data after filtering
      if (nrow(vmr_filtered) > 0) {
        get_enrichment(vmr_filtered, tissue, pop, hfilter)
      } else {
        message("No data left after filtering for ", tissue, " - ", hfilter)
      }
    }
  }
}


# Reproducibility info
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
