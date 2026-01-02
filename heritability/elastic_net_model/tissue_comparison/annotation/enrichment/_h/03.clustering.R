#### Perform hierarchical clustering of enrichment ####

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(data.table)
  library(tidyverse)
  library(dendextend)
})

load_annotation_enrichment <- function(){
  return(data.table::fread("annotation_vmr_enrichment_analysis.txt"))
}

annot <- memoise::memoise(load_annotation_enrichment)

gen_data <- function(){
  err = 0.0000001
  dt <- annot() %>% mutate(across(where(is.character), as.factor)) %>%
    mutate(h2_Category=fct_relevel(h2_Category, rev), `-log10(FDR)`= -log10(FDR),
           `OR Percentile`= OR / (1+OR), p.fdr.sig=FDR < 0.05,
           `log2(OR)` = log2(OR+err),
           p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                         labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                         include.lowest=TRUE))
  return(dt)
}

memDF <- memoise::memoise(gen_data)
df <- memDF() %>% filter(is.finite(`log2(OR)`))

# Exploratory plots
hc <- hclust(dist(df[, 8])) %>%
  as.dendrogram

tissue_cols <- c(
  "Caudate" = "#B36F61",
  "DLPFC" = "#7372A6",
  "Hippocampus" = "#E3C962"
)
labels(hc) <- df$Annotation
labels_colors(hc) <- tissue_cols[df$Tissue]
plot(hc)

cluster <- cutree(hc, k = 3)
table(cluster)

# Brain region dendrogram comparison
dlpfc <- df %>%
  filter(Tissue == "DLPFC")
hc_dlpfc <- as.dendrogram(hclust(dist(dlpfc[, 8])))

h2_cols <- c(
  "Heritable" = "#497C8A",
  "Non-heritable" = "#8CA77B",
  "Low prediction" = "#E3A27F"
)
labels(hc_dlpfc) <- dlpfc$Annotation
labels_colors(hc_dlpfc) <- h2_cols[dlpfc$h2_Category]
plot(hc_dlpfc, cex = 0.7)

caudate <- df %>%
  filter(Tissue == "Caudate")
hc_caudate <- as.dendrogram(hclust(dist(caudate[, 8])))
labels(hc_caudate) <- caudate$Annotation
labels_colors(hc_caudate) <- h2_cols[caudate$h2_Category]
plot(hc_caudate, cex = 0.7)

hippo <- df %>%
  filter(Tissue == "Hippocampus")
hc_hippo <- as.dendrogram(hclust(dist(hippo[, 8])))
labels(hc_hippo) <- hippo$Annotation
labels_colors(hc_hippo) <- h2_cols[hippo$h2_Category]
plot(hc_hippo, cex = 0.7)

dend_diff(hc_caudate, hc_dlpfc)
tanglegram(hc_dlpfc, hc_hippo)

# PCA testing
caudate_pc <- caudate %>%
  select("log2(OR)", "-log10(FDR)")
pc_caudate <- prcomp(caudate_pc, scale = TRUE)
summary(pc_caudate)

plot(pc_caudate$x[, 1], pc_caudate$x[, 2], col = h2_cols[df$h2_Category])

dlpfc_pc <- dlpfc %>%
  select("log2(OR)", "-log10(FDR)")
pc_dlpfc <- prcomp(dlpfc_pc, scale = TRUE)
summary(pc_dlpfc)

plot(pc_dlpfc$x[, 1], pc_dlpfc$x[, 2], col = cols[df$h2_Category])

hippo_pc <- hippo %>%
  select("log2(OR)", "-log10(FDR)")
pc_hippo <- prcomp(hippo_pc, scale = TRUE)
summary(pc_hippo)

plot(pc_hippo$x[, 1], pc_hippo$x[, 2], col = cols[df$h2_Category])

caudate_pc_test <- caudate_vmrs %>%
  select_if(is.numeric) %>%
  select(-start, -end)

pc_caudate <- prcomp(caudate_pc_test)
summary(pc_caudate)
plot(pc_caudate$x[, 1], pc_caudate$x[, 2])

# UMAP testing

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()