## This script maps shared VMRs for mash modeling input

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(tidyr)
  library(dplyr)
})

## Main
# Create output dir
output_path <- here("meqtl-analysis", "vmrs", "shared_vmrs", "_m")

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# Get fractional threshold from command line
args <- commandArgs(trailingOnly = TRUE)
option_flag <- args[1]

# Read in bed files
caudate_bed_fn <- here("meqtl-analysis", "vmrs", "caudate", "_m", "feature.bed")
hippo_bed_fn <- here("meqtl-analysis", "vmrs", "hippocampus", "_m", "feature.bed")
dlpfc_bed_fn <- here("meqtl-analysis", "vmrs", "dlpfc", "_m", "feature.bed")
caudate_hippo_fn <- here("heritability", "elastic_net_model", "tissue_comparison",
                         "upset_plot", "_m", option_flag, "sets", "caudate_hippocampus_overlap_all.bed")
caudate_dlpfc_fn <- here("heritability", "elastic_net_model", "tissue_comparison",
                      "upset_plot", "_m", option_flag, "sets", "caudate_dlpfc_overlap_all.bed")

# Add binary values for each tissue
caudate_bed <- fread(caudate_bed_fn, select=2:5,
                     col.names = c("seqnames_caudate", "start_caudate", 
                                   "end_caudate", "feature_id_caudate")) %>%
  mutate(caudate = 1L)
hippo_bed <- fread(hippo_bed_fn, select=2:5,
                     col.names = c("seqnames_hippocampus", "start_hippocampus", 
                                   "end_hippocampus", "feature_id_hippocampus"))
dlpfc_bed <- fread(dlpfc_bed_fn, select=2:5,
                     col.names = c("seqnames_dlpfc", "start_dlpfc", 
                                   "end_dlpfc", "feature_id_dlpfc"))

# Group shared VMRs
caudate_hippo <- fread(caudate_hippo_fn, select=1:6, 
                       col.names = c("seqnames_caudate", "start_caudate", "end_caudate",
                                     "seqnames_hippocampus", "start_hippocampus", "end_hippocampus")) %>%
  mutate(seqnames_caudate = paste0("chr", seqnames_caudate), 
         seqnames_hippocampus = paste0("chr", seqnames_hippocampus), hippocampus = 1L) %>%
  left_join(hippo_bed, by = c("seqnames_hippocampus", "start_hippocampus", "end_hippocampus")) # Add region-specific feature IDs
caudate_dlpfc <- fread(caudate_dlpfc_fn, select=1:6,
                       col.names = c("seqnames_caudate", "start_caudate", "end_caudate",
                                     "seqnames_dlpfc", "start_dlpfc", "end_dlpfc")) %>%
  mutate(seqnames_caudate = paste0("chr", seqnames_caudate), 
         seqnames_dlpfc = paste0("chr", seqnames_dlpfc), dlpfc = 1L) %>%
  left_join(dlpfc_bed, by = c("seqnames_dlpfc", "start_dlpfc", "end_dlpfc")) # Add region-specific feature IDs

# Merge vmr overlaps
merged <- caudate_bed %>%
  left_join(caudate_hippo, by = c("seqnames_caudate", "start_caudate", "end_caudate")) %>%
  left_join(caudate_dlpfc, by = c("seqnames_caudate", "start_caudate", "end_caudate")) %>%
  mutate(hippocampus = ifelse(is.na(hippocampus), 0L, hippocampus),
         dlpfc = ifelse(is.na(dlpfc), 0L, dlpfc)) 

# Add shared key for each cluster
shared <- merged %>%
  filter(caudate == 1 & hippocampus == 1 & dlpfc == 1) %>%
  mutate(shared_feature_id = paste0("VMR_S", row_number())) %>%
  select(shared_feature_id, everything())

# Write for mash input
out_mash <- file.path(output_path, "TOPMed_LIBD_shared_vmr_key.tsv")
fwrite(shared, out_mash, sep='\t', row.names = FALSE)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()