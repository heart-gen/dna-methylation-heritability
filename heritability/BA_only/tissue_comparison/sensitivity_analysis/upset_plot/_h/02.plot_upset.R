## This script generates the upset plot for VMRs across brain regions

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(UpSetR)
  library(ggplot2)
  #library(ComplexUpset)
  library(tibble)
})

## Function
count_intersections <- function(fn){
    vmrs <- fread(fn)
    return(nrow(vmrs))
}

# Formatting

# Plot Upset
plot_upset <- function(sets, fn){
  pdf(fn, width = 6, height = 4)
  upset <- UpSetR::upset(
    fromExpression(sets),
    order.by = "freq",
    sets = c("Caudate", "DLPFC", "Hippocampus")
  )
  print(upset)
  dev.off()
}

# Testing complex upset
#df <- as.data.frame(fromExpression(sets))
#ComplexUpset::upset(df, intersect = c("Caudate", "DLPFC", "Hippocampus"))

## Main
out_path <- here("heritability/elastic_net_model/tissue_comparison/upset_plot/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

option_flags <- c("F_0.25", "F_0.5", "F_0.75", "f_0.25", "f_0.5", "f_0.75")
h2_categories <- c("heritable", "non-heritable", "low_prediction", "all")

for (flag in option_flags) {
  for (h2_cat in h2_categories){
    
    # Get overlapping VMRs
    sets <- c(
      Caudate = count_intersections(paste0("./", flag, "/sets/caudate_specific_", h2_cat, ".bed")),
      DLPFC = count_intersections(paste0("./", flag, "/sets/dlpfc_specific_", h2_cat, ".bed")),
      Hippocampus = count_intersections(paste0("./", flag, "/sets/hippocampus_specific_", h2_cat, ".bed")),
      `Caudate&Hippocampus` = count_intersections(paste0("./", flag, "/sets/caudate_hippocampus_overlap_", h2_cat, ".bed")),
      `Caudate&DLPFC` = count_intersections(paste0("./", flag, "/sets/caudate_dlpfc_overlap_", h2_cat, ".bed")),
      `Hippocampus&DLPFC` = count_intersections(paste0("./", flag, "/sets/hippocampus_dlpfc_overlap_", h2_cat, ".bed")),
      `Caudate&Hippocampus&DLPFC` = count_intersections(paste0("./", flag, "/sets/3tissues_overlap_", h2_cat, ".bed.tmp"))
    )
    
    # Plot upset
    out_upset <- file.path(out_path, flag, paste0("VMR_upsetR_", h2_cat, ".pdf"))
    plot_upset(sets, out_upset)

  }
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
