#### Perform hierarchical clustering of enrichment ####

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(data.table)
  library(tidyverse)
  library(dendextend)
})

# Function 
load_annotation_enrichment <- function(){
  return(data.table::fread("annotation_vmr_enrichment_analysis.txt"))
}

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

plot_dendrogram <- function(df, tissue, h2_cols, output_path) {
  tissue_df <- df %>%
    filter(Tissue == tissue)
  hc <- as.dendrogram(hclust(dist(tissue_df[, 8])))
  
  h2_cols <- c(
    "Heritable" = "#497C8A",
    "Non-heritable" = "#8CA77B",
    "Low prediction" = "#E3A27F"
  )
  labels(hc) <- tissue_df$Annotation
  labels_colors(hc) <- h2_cols[tissue_df$h2_Category]
  
  pdf(file = file.path(output_path, paste0(tissue, "_hc.pdf")))
  plot(hc, cex = 0.7)
  dev.off()
  
  return(hc)
}

pca <- function(df, tissue, h2_cols, output_path) {
  filtered <- df %>%
    filter(Tissue == tissue) %>%
    select("log2(OR)", "-log10(FDR)")
  pc <- prcomp(filtered, scale = TRUE)
  print(summary(pc))
  
  pdf(file = file.path(output_path, paste0(tissue, "_pca.pdf")))
  plot(pc$x[, 1], pc$x[, 2], col = h2_cols[df$h2_Category])
  dev.off()
  
  return(pc)
}

# Main
annot <- memoise::memoise(load_annotation_enrichment)
memDF <- memoise::memoise(gen_data)
df <- memDF() %>% filter(is.finite(`log2(OR)`))

output_path <- here("heritability", "elastic_net_model", "tissue_comparison",
                    "annotation", "enrichment", "_m")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

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
tissues <- c("Caudate", "DLPFC", "Hippocampus")
h2_cols <- c(
  "Heritable" = "#497C8A",
  "Non-heritable" = "#8CA77B",
  "Low prediction" = "#E3A27F"
)

for (tissue in tissues) {
  hc <- plot_dendrogram(df, tissue, h2_cols, output_path)
  pc <- pca(df, tissue, h2_cols, output_path)
}

dend_diff(hc_caudate, hc_dlpfc)
tanglegram(hc_dlpfc, hc_hippo)

# UMAP testing

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()