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
dlpfc_bed_fn <- here("meqtl-analysis", "vmrs", "dlpfc", "_m", "feature.bed")
dlpfc_hippo_fn <- here("heritability", "elastic_net_model", "tissue_comparison",
                         "upset_plot", "_m", "f_0.25", "sets", "hippocampus_dlpfc_overlap_all.bed")
dlpfc_caudate_fn <- here("heritability", "elastic_net_model", "tissue_comparison",
                      "upset_plot", "_m", "f_0.25", "sets", "caudate_dlpfc_overlap_all.bed")

# Add binary values for each tissue
dlpfc_bed <- fread(dlpfc_bed_fn, select=2:5) %>%
  mutate(caudate = 1L)
dlpfc_hippo <- fread(dlpfc_hippo_fn, select=4:6, 
                       col.names = c("seqnames", "start", "end")) %>%
  mutate(seqnames = paste0("chr", seqnames), hippocampus = 1L)
dlpfc_caudate <- fread(dlpfc_caudate_fn, select=4:6,
                       col.names = c("seqnames", "start", "end")) %>%
  mutate(seqnames = paste0("chr", seqnames), caudate = 1L)

# Merge vmr overlaps
merged <- dlpfc_bed %>%
  left_join(dlpfc_hippo, by = c("seqnames", "start", "end")) %>%
  left_join(dlpfc_caudate, by = c("seqnames", "start", "end")) %>%
  mutate(hippocampus = ifelse(is.na(hippocampus), 0L, hippocampus),
         caudate = ifelse(is.na(caudate), 0L, caudate)) 

# Add shared key for each cluster
merged <- merged %>%
  mutate(
    vmr_key = case_when(
      hippocampus == 0 & caudate == 0 ~ "dlpfc",
      hippocampus == 1 & caudate == 0 ~ "dlpfc_hippocampus",
      hippocampus == 0 & caudate == 1 ~ "dlpfc_caudate",
      hippocampus == 1 & caudate == 1 ~ "caudate_hippocampus_dlpfc"
    )
  )

# Write for mash input
out_mash <- file.path(output_path, "shared_vmr_key.tsv")
fwrite(merged, out_mash, sep='\t', row.names = FALSE)
