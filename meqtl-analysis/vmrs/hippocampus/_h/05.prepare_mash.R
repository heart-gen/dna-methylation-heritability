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
hippo_bed_fn <- here("meqtl-analysis", "vmrs", "hippocampus", "_m", "feature.bed")
hippo_caudate_fn <- here("heritability", "elastic_net_model", "tissue_comparison",
                         "upset_plot", "_m", "f_0.25", "sets", "caudate_hippocampus_overlap_all.bed")
hippo_dlpfc_fn <- here("heritability", "elastic_net_model", "tissue_comparison",
                      "upset_plot", "_m", "f_0.25", "sets", "hippo_dlpfc_overlap_all.bed")

# Add binary values for each tissue
hippo_bed <- fread(hippo_bed_fn, select=2:5) %>%
  mutate(hippocampus = 1L)
hippo_caudate <- fread(hippo_caudate_fn, select=4:6, 
                       col.names = c("seqnames", "start", "end")) %>%
  mutate(seqnames = paste0("chr", seqnames), caudate = 1L)
hippo_dlpfc <- fread(hippo_dlpfc_fn, select=1:3,
                       col.names = c("seqnames", "start", "end")) %>%
  mutate(seqnames = paste0("chr", seqnames), dlpfc = 1L)

# Merge vmr overlaps
merged <- hippo_bed %>%
  left_join(hippo_caudate, by = c("seqnames", "start", "end")) %>%
  left_join(hippo_dlpfc, by = c("seqnames", "start", "end")) %>%
  mutate(caudate = ifelse(is.na(caudate), 0L, caudate),
         dlpfc = ifelse(is.na(dlpfc), 0L, dlpfc)) 

# Add shared key for each cluster
merged <- merged %>%
  mutate(
    vmr_key = case_when(
      caudate == 0 & dlpfc == 0 ~ "hippocampus",
      caudate == 1 & dlpfc == 0 ~ "hippocampus_caudate",
      caudate == 0 & dlpfc == 1 ~ "hippocampus_dlpfc",
      caudate == 1 & dlpfc == 1 ~ "hippocampus_hippocampus_dlpfc"
    )
  )

# Write for mash input
out_mash <- file.path(output_path, "shared_vmr_key.tsv")
fwrite(merged, out_mash, sep='\t', row.names = FALSE)
