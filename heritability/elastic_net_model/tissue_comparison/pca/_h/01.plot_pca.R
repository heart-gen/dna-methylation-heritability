#### Plot PCA of VMRs ####

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(data.table)
  library(tidyverse)
  library(ggplot2)
})

# Function 
get_vmrs <- function(tissue){
  enet_file <- here("heritability/elastic_net_model/", 
                    paste0(tolower(tissue), "/_m/", tolower(tissue), "_summary_elastic-net.tsv"))
  vmr <- data.table::fread(enet_file) %>%
    mutate(chrom = paste0("chr", chrom)) %>%
    rename(chr = chrom) %>%
    na.omit() %>%
    mutate(h2_category = case_when(
      r_squared_cv <= 0.75 ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
    ))
  
  return(vmr)
}

save_plot <- function(p, fn, w, h, dpi){
  for(ext in c(".pdf", ".png")){
    ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h, dpi=dpi)
  }
}

pca <- function(df, tissue, h2_cols, output_path) {
  filtered <- df %>%
    filter(region == tissue) %>%
    select("h2_unscaled", "r_squared_cv", h2_category)
  pc <- prcomp(filtered[, c("h2_unscaled", "r_squared_cv")], scale = TRUE)
  print(summary(pc))

  # format df for plotting 
  pc_df <- as.data.frame(pc$x[ , 1:2])
  pc_df <- cbind(pc_df, filtered$h2_category)
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
output_path <- here("heritability", "elastic_net_model", "tissue_comparison",
                    "pca", "_m")
if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

tissue_cols <- c(
  "Caudate" = "#B36F61",
  "DLPFC" = "#7372A6",
  "Hippocampus" = "#E3C962"
)

tissues <- c("caudate", "DLPFC", "Hippocampus")
h2_cols <- c(
  "Heritable" = "#497C8A",
  "Non-heritable" = "#8CA77B",
  "Low prediction" = "#E3A27F"
)

for (tissue in tissues) {
  vmr <- get_vmrs(tissue)
  p <- pca(vmr, tissue, h2_cols, output_path)
  pca_fn <- file.path(output_path, paste0(tolower(tissue), "_pca"))
  save_plot(p, pca_fn, w = 10, h = 6, dpi = 300)
}

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()