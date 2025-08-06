#### Get genes nearby VMRs ####
suppressPackageStartupMessages({
    library(dplyr)
    library(GenomicRanges)
    library(here)
    library(data.table)
    library(bumphunter)
    library(rtracklayer)
    library(txdbmaker)
})

# Function
load_vmrs <- function(tissue) {
                                        # combine vmr bed files from all chr
  vmr_files <- here(paste0("heritability/", tissue, "/_m/vmr/chr_", 1:22, 
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

load_annotations <- function(annotation_fn) {
  
                                        # load gene annotations 
  txdb <- makeTxDbFromGFF(annotation_fn, format = "gtf")
  annotated_genes <- annotateTranscripts(txdb, annotationPackage = "org.Hs.eg.db")
  
  return(annotated_genes)
}

annotate_vmrs <- function(vmr_gr, annotated_genes, out_file) {
  
                                        # find nearest genes/transcripts
  matched <- matchGenes(vmr_gr, annotated_genes)
  
                                        # extract vmr positions                                          
  vmr     <- data.table(chr = as.character(seqnames(vmr_gr)),
                        start = start(vmr_gr),
                        end = end(vmr_gr))
  
                                        # merge annotations with vmr positions  
  matched <- cbind(vmr, matched)
  
                                        # write annotations to file
  fwrite(matched, out_file, sep = "\t")
  return(matched)
}

# Main
  
                                        # create output dir if it doesn't exist
out_path <- here("heritability/elastic_net_model/tissue_comparison/annotation/_m")

if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

annotation_fn <- "/projects/b1213/resources/genomes/human/gencode-v47/gtf/gencode.v47.primary_assembly.annotation.gtf"

tissues <- c("dlpfc", "caudate", "hippocampus")

for (tissue in tissues) {
  out_file         <- file.path(out_path, 
                                paste0(tissue, "_vmr_gene_annotation_match_genes.tsv"))
  vmr_gr           <- load_vmrs(tissue)
  annotated_genes  <- load_annotations(annotation_fn)
  matched          <- annotate_vmrs(vmr_gr, annotated_genes, out_file)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()