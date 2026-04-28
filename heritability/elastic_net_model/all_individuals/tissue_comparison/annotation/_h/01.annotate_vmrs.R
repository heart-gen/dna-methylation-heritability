#### Get genes nearby VMRs ####
suppressPackageStartupMessages({
    library(dplyr)
    library(GenomicRanges)
    library(here)
    library(data.table)
    library(bumphunter)
    library(rtracklayer)
    library(txdbmaker)
    library(annotatr)
    library(tidyr)
})

## Function
load_vmrs <- function(tissue) {
                                        # combine vmr bed files from all chr
  vmr_files <- here(paste0("vmr-analysis/all_individuals/", tissue, 
                    "/_m/vmr/chr_", c(1:22, "X", "Y"), "/vmr.bed"))
  
  vmr_list  <- lapply(vmr_files, function(f) {
    if (file.exists(f)) {
      fread(f, col.names = c("chr", "start", "end"))
    } else {
      message("File does not exist: ", f)
      return(NULL)
    }
  })
  vmr_list <- Filter(Negate(is.null), vmr_list)
  vmr      <- rbindlist(vmr_list)
  
                                        # convert to GRanges
  vmr$chr  <- ifelse(grepl("^chr", vmr$chr), vmr$chr, paste0("chr", vmr$chr))
  vmr_gr   <- GRanges(seqnames = vmr$chr,
                      ranges = IRanges(start = vmr$start, end = vmr$end))
  
  return(vmr_gr)
  
}

load_annotations <- function() {
  
                                        # load genic annotations  
  genic <- build_annotations(
    genome = "hg38",
    annotations = "hg38_basicgenes"
  )
  
                                        # load intergenic annotations
  intergenic <- build_annotations(
    genome = "hg38",
    annotations = "hg38_genes_intergenic"
  )
  
                                        # load cpg annotations
  cpg_islands <- build_annotations(
    genome = "hg38",
    annotations = "hg38_cpgs"
  )
  
                                        # load FANTOM5 enhancer annotations
  enhancers <- build_annotations(
    genome = "hg38",
    annotations = "hg38_enhancers_fantom"
  )
  
                                        # compile annotations
  annots <- c(genic, intergenic, cpg_islands, enhancers)
  
  return(annots)
}

annotate_vmrs <- function(vmr_gr, enet, annots, out_file) {
  
                                        # annotate vmrs
  annotated_vmrs <- annotate_regions(
    regions = vmr_gr,
    annotations = annots,
    ignore.strand = TRUE)
  
                                        # get h2 categories
  enet <- na.omit(enet)
  enet <- enet %>% 
    dplyr::select(chrom, start, end, h2_category) %>% 
    dplyr::rename("chr" = "chrom") %>%
    mutate(chr = paste0("chr", chr)) %>%
    h2_category = factor(h2_category,
                         levels = c("Heritable", "Non-heritable", "Low prediction"))
  
                                        # merge annotations with h2  
  annot_df <- data.frame(annotated_vmrs)
  merged   <- annot_df %>%
    left_join(enet, by = c("seqnames" = "chr", "start", "end"))
  
                                        # write annotations to file
  fwrite(merged, out_file, sep = "\t")
  return(merged)
}

format_annotations <- function(merged, out_file) {
  
                                        # pivot wide
  wide <- merged %>%
    dplyr::select(seqnames, start, end, annot.type, h2_category) %>%
    mutate(val = 1) %>%
    pivot_wider(
      names_from = annot.type,
      values_from = val, 
      values_fn = max,
      values_fill = 0     
    )
                                        # write formatted annotations
  fwrite(wide, out_file, sep = "\t")
}

## Main
  
                                        # create output dir if it doesn't exist
out_path <- here("heritability/elastic_net_model/all_individuals/tissue_comparison/annotation/_m")

if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

tissues <- c("dlpfc", "caudate", "hippocampus")

for (tissue in tissues) {
  out_annot  <- file.path(out_path, 
                         paste0(tissue, "_vmr_annotations_hg38.tsv"))
  out_wide   <- file.path(out_path, 
                          paste0(tissue, "_vmr_annotations_hg38_wide.tsv"))
  vmr_gr    <- load_vmrs(tissue)
  annots    <- load_annotations()
  
  enet_file <- here("heritability/elastic_net_model/all_individuals", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net_matched_r2_0.3.tsv"))
  enet      <- read.table(enet_file, sep = "\t", header = TRUE)
  
  merged    <- annotate_vmrs(vmr_gr, enet, annots, out_annot)
  
  summary <- merged %>%
    group_by(h2_category, annot.type) %>%
    summarise(n = n())
  print(summary)
  
  wide <- format_annotations(merged, out_wide)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
