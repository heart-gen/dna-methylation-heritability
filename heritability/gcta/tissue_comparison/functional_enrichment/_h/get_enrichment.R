##### Performs GO and KEGG enrichment analysis for all brain regions #####
suppressPackageStartupMessages({
    library(rGREAT) 
    library(here)
    library(dplyr)
})

## Function
get_enrichment <- function(tissue, enrichment_type, heritability_filter) {
    # Read in summary table 
    vmr <- read.table(here("heritability/gcta/", paste0(region, "/_m/summary/greml_summary.tsv"), sep="\t", header=TRUE))

    # Filter by heritability
    if (heritability_filter == "heritable") {
        vmr_filtered <- vmr %>% filter(Sum.of.V.G._Vp_Variance > 0.05)
    } else if (heritability_filter == "non_heritable") {
        vmr_filtered <- vmr %>% filter(Sum.of.V.G._Vp_Variance <= 0.05)
    } else {
        vmr_filtered <- vmr
    }

    colnames(vmr_filtered) <- c("seqnames","start","end")
    vmr <- plyranges::as_granges(vmr_filtered)
   
    if (enrichment_type == "go") {
        # GO enrich
        res <- great(vmr, "GO:BP", "RefSeq:hg38", background=as.character(1:22))
        tb$GO  <- getEnrichmentTable(res)
    } else if (enrichment_type == "KEGG") {
        # KEGG enrich
        res <- great(vmr, "kegg", "RefSeq:hg38", background=as.character(1:22))
        tb$KEGG <- getEnrichmentTable(res)
    }

    write.csv(tb[[enrichment_type]], paste0(region, "_", heritability_filter, "_", enrichment_type, ".csv"), row.names = FALSE)
}

## Main
enrichment <- c("go", "kegg")
tissues <- c("caudate", "dlpfc", "hippocampus")
heritability <- c("all", "heritable", "non_heritable")

for (tissue in tissues) {
  for (enrichment_type in enrichment) {
    for (heritability_filter in heritability) {
      get_enrichment(tissue, enrichment_type, heritability_filter)
    }
  }
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()