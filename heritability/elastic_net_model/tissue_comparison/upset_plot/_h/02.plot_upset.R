## This script generates the upset plot for VMRs across brain regions

suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(ComplexHeatmap)
  library(UpSetR)
  library(ggplot2)
  library(ComplexUpset)
  library(tibble)
})

## Function
count_intersections <- function(fn){
    vmrs <- fread(fn)
    return(nrow(vmrs))
}

# Formatting

# Plot Upset
plot_upset <- function(sets){
  pdf("VMR_upsetR_plot.pdf", width = 6, height = 4)
  UpSetR::upset(
    fromExpression(sets),
    order.by = "freq",
    sets = c("Caudate", "DLPFC", "Hippocampus")
  )
  dev.off()
}

# Testing complex upset
df <- as.data.frame(fromExpression(sets))
ComplexUpset::upset(df, intersect = c("Caudate", "DLPFC", "Hippocampus"))

## Main
sets <- c(
  Caudate = count_intersections("./f_0.25/sets/caudate_specific.bed"),
  DLPFC = count_intersections("./f_0.25/sets/dlpfc_specific.bed"),
  Hippocampus = count_intersections("./f_0.25/sets/hippocampus_specific.bed"),
  `Caudate&Hippocampus` = count_intersections("./f_0.25/sets/caudate_hippocampus_overlap.bed"),
  `Caudate&DLPFC` = count_intersections("./f_0.25/sets/caudate_dlpfc_overlap.bed"),
  `Hippocampus&DLPFC` = count_intersections("./f_0.25/sets/hippocampus_dlpfc_overlap.bed"),
  `Caudate&Hippocampus&DLPFC` = count_intersections("./f_0.25/sets/3tissues_overlap.bed.tmp")
)

plot_upset(sets)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
