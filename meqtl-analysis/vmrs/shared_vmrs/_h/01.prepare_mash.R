## This script maps shared VMRs for mash modeling input

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(tidyr)
  library(dplyr)
})

## Main
get_bed <- function(tissue) {
  bed_fn <- here("meqtl-analysis", "vmrs", tissue, "_m", "feature.bed")
  bed <- fread(bed_fn, select=2:5,
               col.names = c(paste0("seqnames_", tissue), paste0("start_", tissue), 
                             paste0("end_", tissue), paste0("feature_id_", tissue)))
  return(bed)
}

get_overlap <- function(tissue1, tissue2, option_flag){
  overlap_bed_fn <- here("heritability", "elastic_net_model", "tissue_comparison",
                         "upset_plot", "_m", option_flag, "sets", 
                         paste0(tissue1, "_", tissue2, "_overlap_all.bed"))
  overlap_bed <- fread(overlap_bed_fn, select=1:6,
               col.names = c(paste0("seqnames_", tissue1), paste0("start_", tissue1), 
                             paste0("end_", tissue1), paste0("seqnames_", tissue2), 
                             paste0("start_", tissue2), paste0("end_", tissue2)))
  return(overlap_bed)
}

write_shared_vmr <- function(merged, tissues, output_path){
  # Add shared key for each cluster
  shared <- merged %>%
    filter(if_all(all_of(tissues), ~ .x == 1L)) %>%
    mutate(shared_feature_id = paste0("VMR_S", row_number())) %>%
    select(shared_feature_id, everything())
  
  # Write for mash input
  out_mash <- file.path(output_path, "TOPMed_LIBD_shared_vmr_key.tsv")
  fwrite(shared, out_mash, sep='\t', row.names = FALSE)
  
  return(shared)
}

## Main
# Create output dir
output_path <- here("meqtl-analysis", "vmrs", "shared_vmrs", "_m")

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# Get fractional threshold from command line
args <- commandArgs(trailingOnly = TRUE)
option_flag <- args[1]
tissues <- c("caudate", "hippocampus", "dlpfc")

# Read in bed files
caudate_bed <- get_bed("caudate") %>%
  mutate(caudate = 1L) #anchor tissue
hippo_bed <- get_bed("hippocampus")
dlpfc_bed <- get_bed("dlpfc")

# Group shared VMRs
caudate_hippo <- get_overlap("caudate", "hippocampus", option_flag) %>%
  mutate(seqnames_caudate = paste0("chr", seqnames_caudate), 
         seqnames_hippocampus = paste0("chr", seqnames_hippocampus), 
         hippocampus = 1L) %>%
  left_join(hippo_bed, 
            by = c("seqnames_hippocampus", "start_hippocampus", "end_hippocampus")) # Add region-specific feature IDs
caudate_dlpfc <- get_overlap("caudate", "dlpfc", option_flag) %>%
  mutate(seqnames_caudate = paste0("chr", seqnames_caudate), 
         seqnames_dlpfc = paste0("chr", seqnames_dlpfc), 
         dlpfc = 1L) %>%
  left_join(dlpfc_bed, 
            by = c("seqnames_dlpfc", "start_dlpfc", "end_dlpfc")) # Add region-specific feature IDs

# Merge vmr overlaps
merged <- caudate_bed %>%
  left_join(caudate_hippo, by = c("seqnames_caudate", "start_caudate", "end_caudate")) %>%
  left_join(caudate_dlpfc, by = c("seqnames_caudate", "start_caudate", "end_caudate")) %>%
  mutate(hippocampus = ifelse(is.na(hippocampus), 0L, hippocampus),
         dlpfc = ifelse(is.na(dlpfc), 0L, dlpfc)) 

# Write to file
shared <- write_shared_vmr(merged, tissues, output_path)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()