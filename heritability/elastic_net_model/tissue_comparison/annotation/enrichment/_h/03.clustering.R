#### Perform hierarchical clustering of enrichment ####

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(data.table)
  library(tidyverse)
  library(dendextend)
  library(M3C)
  library(ggplot2)
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

save_plot <- function(p, fn, w, h, dpi){
  for(ext in c(".pdf", ".png")){
    ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h, dpi=dpi)
  }
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
  
  pdf(file = file.path(output_path, "hc", paste0(tolower(tissue), "_hc.pdf")))
  plot(hc, cex = 0.7)
  dev.off()
  
  return(hc)
}

pca <- function(df, tissue, h2_cols, output_path) {
  filtered <- df %>%
    filter(Tissue == tissue) %>%
    select("log2(OR)", "-log10(FDR)", h2_Category)
  pc <- prcomp(filtered[, c("log2(OR)", "-log10(FDR)")], scale = TRUE)
  print(summary(pc))

  # format df for plotting 
  pc_df <- as.data.frame(pc$x[ , 1:2])
  pc_df <- cbind(pc_df, filtered$h2_Category)
  colnames(pc_df) <- c("PC1", "PC2", "Heritability Category")
  
  # plot pca
  p <- pc_df %>%
    ggplot(aes(x = PC1, y = PC2, 
               color = `Heritability Category`, 
               fill = `Heritability Category`, 
               shape = `Heritability Category`)) +
    geom_point(size = 5) +
    labs(x = "PC1",
         y = "PC2",
         color = "Heritability Category",
         shape = "Heritability Category") +
    scale_fill_manual(values = h2_cols) +
    scale_color_manual(values = h2_cols) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 1, linetype = "solid")
    )
  
  return(p)
}

# Main
annot <- memoise::memoise(load_annotation_enrichment)
memDF <- memoise::memoise(gen_data)
df <- memDF() %>% filter(is.finite(`log2(OR)`))

output_path <- here("heritability", "elastic_net_model", "tissue_comparison",
                    "annotation", "enrichment", "_m")
subdirs <- c("hc", "pca")
for (subdir in subdirs){
  subdir_path <- file.path(output_path, subdir)
  if (!dir.exists(subdir_path)) {
    dir.create(subdir_path, recursive = TRUE)
  }
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
  p <- pca(df, tissue, h2_cols, output_path)
  pca_fn <- file.path(output_path, "pca", paste0(tolower(tissue), "_pca"))
  save_plot(p, pca_fn, w = 10, h = 6, dpi = 300)
}

#dend_diff(hc_caudate, hc_dlpfc)
#tanglegram(hc_dlpfc, hc_hippo)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()