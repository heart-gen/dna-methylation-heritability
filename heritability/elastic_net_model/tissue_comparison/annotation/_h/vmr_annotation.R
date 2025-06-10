#### Get genes nearby VMRs ####
suppressPackageStartupMessages({
    library(dplyr)
    library(GenomicRanges)
    library(biomaRt)
    library(here)
    library(data.table)
})

# Function
load_vmrs <- function(tissue) {
                                        # combine vmr bed files from all chr
  vmr_files <- here(paste0("heritability/gcta/", tissue, "/_m/vmr/chr_", 1:22, 
                           "/vmr.bed"))
  
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
                                        # load gene annotations 
  ensembl <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")
  genes   <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 
                                'start_position', 'end_position'), 
                   mart=ensembl)
  genes   <- genes[genes$chromosome_name %in% c(1:22), ]
  
                                        # convert to GRanges
  gene_gr <- GRanges(seqnames  = paste0("chr", genes$chromosome_name),
                     ranges    = IRanges(genes$start_position, genes$end_position),
                     gene_id   = genes$ensembl_gene_id,
                     gene_name = genes$external_gene_name)
  
  return(gene_gr)
}

annotate_vmrs <- function(vmr_gr, gene_gr, out_file) {
  
                                        # find nearest genes 
  nearest_genes   <- distanceToNearest(vmr_gr, gene_gr)
  nearest_gene_id <- gene_gr$gene_id[subjectHits(nearest_genes)]
  distance        <- mcols(nearest_genes)$distance
                                          
  vmr             <- data.table(chr = as.character(seqnames(vmr_gr)),
                                start = start(vmr_gr),
                                end = end(vmr_gr),
                                nearest_gene_id = nearest_gene_id,
                                distance = distance)
  
                                        # write annotations to file
  fwrite(vmr, out_file, sep = "\t")
  return(vmr)
}

# Main
  
                                        # create output dir if it doesn't exist
out_path <- here("heritability/elastic_net_model/tissue_comparison/annotation/_m")

if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

tissues <- c("dlpfc", "caudate", "hippocampus")

for (tissue in tissues) {
  out_file <- file.path(out_path, paste0(tissue, "_vmr_gene_annotation.tsv"))
  vmr_gr   <- load_vmrs(tissue)
  gene_gr  <- load_annotations()
  vmr      <- annotate_vmrs(vmr_gr, gene_gr, out_file)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()