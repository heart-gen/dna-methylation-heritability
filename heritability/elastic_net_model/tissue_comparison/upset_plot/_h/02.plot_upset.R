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
plot_upset <- function(sets, fn, w, h){
  for(ext in c('.pdf')){
    pdf(file=paste0(fn, ext), width=w, height=h)
    grid::grid.newpage()
    UpSetR::upset(
      fromExpression(sets),
      order.by = "freq",
      sets = c("Caudate", "DLPFC", "Hippocampus")
    )
    dev.off()
  }
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

for (flag in option_flags) {
  
  # Get overlapping VMRs
  sets <- c(
    Caudate = count_intersections(paste0("./", flag, "/sets/caudate_specific.bed")),
    DLPFC = count_intersections(paste0("./", flag, "/sets/dlpfc_specific.bed")),
    Hippocampus = count_intersections(paste0("./", flag, "/sets/hippocampus_specific.bed")),
    `Caudate&Hippocampus` = count_intersections(paste0("./", flag, "/sets/caudate_hippocampus_overlap.bed")),
    `Caudate&DLPFC` = count_intersections(paste0("./", flag, "/sets/caudate_dlpfc_overlap.bed")),
    `Hippocampus&DLPFC` = count_intersections(paste0("./", flag, "/sets/hippocampus_dlpfc_overlap.bed")),
    `Caudate&Hippocampus&DLPFC` = count_intersections(paste0("./", flag, "/sets/3tissues_overlap.bed.tmp"))
  )
  
  # Plot upset
  out_upset <- file.path(out_path, flag, "VMR_upsetR")
  plot_upset(sets, out_upset, 6, 4)

}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
