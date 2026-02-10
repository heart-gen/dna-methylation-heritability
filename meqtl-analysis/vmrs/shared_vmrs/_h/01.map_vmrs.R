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

option_flags <- c("F_0.25")
h2_categories <- c("all")

# Read in bed files
caudate_bed_fn <- here("meqtl-analysis", "vmrs", "caudate", "_m", "feature.bed")
caudate_hippo_fn <- here("heritability", "elastic_net_model", "tissue_comparison",
                         "upset_plot", "_m", "f_0.25", "sets", "caudate_hippocampus_overlap_all.bed")
caudate_dlpfc_fn <- here("heritability", "elastic_net_model", "tissue_comparison",
                      "upset_plot", "_m", "f_0.25", "sets", "caudate_dlpfc_overlap_all.bed")

# Add binary values for each tissue
caudate_bed <- fread(caudate_bed_fn, select=2:5) %>%
  mutate(caudate = 1L)
caudate_hippo <- fread(caudate_hippo_fn, select=1:3, 
                       col.names = c("seqnames", "start", "end")) %>%
  mutate(seqnames = paste0("chr", seqnames), hippocampus = 1L)
caudate_dlpfc <- fread(caudate_dlpfc_fn, select=1:3,
                       col.names = c("seqnames", "start", "end")) %>%
  mutate(seqnames = paste0("chr", seqnames), dlpfc = 1L)

# Merge vmr overlaps
merged <- caudate_bed %>%
  left_join(caudate_hippo, by = c("seqnames", "start", "end")) %>%
  left_join(caudate_dlpfc, by = c("seqnames", "start", "end")) %>%
  mutate(hippocampus = ifelse(is.na(hippocampus), 0L, hippocampus),
         dlpfc = ifelse(is.na(dlpfc), 0L, dlpfc)) 

# Add shared key for each cluster
merged <- merged %>%
  mutate(
    vmr_key = case_when(
      hippocampus == 0 & dlpfc == 0 ~ "caudate",
      hippocampus == 1 & dlpfc == 0 ~ "caudate_hippocampus",
      hippocampus == 0 & dlpfc == 1 ~ "caudate_dlpfc",
      hippocampus == 1 & dlpfc == 1 ~ "caudate_hippocampus_dlpfc"
    )
  )

# Write for mash input
out_mash <- file.path(output_path, "shared_vmr_key.tsv")
fwrite(merged, out_mash, sep='\t', row.names = FALSE)