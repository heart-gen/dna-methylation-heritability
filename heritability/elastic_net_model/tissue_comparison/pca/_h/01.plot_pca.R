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
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
  vmr <- data.table::fread(enet_file) %>%
    na.omit() %>%
    mutate(h2_category = case_when(
      r_squared_cv <= 0.75 ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
    ))
  
  return(vmr)
}

merge_meth <- function(meth_files){
  meth_list <- vector("list", length(meth_files))
  for (i in seq_along(meth_files)) {
    file  <- meth_files[i]
    # Get pos from filename
    chr   <- basename(dirname(file))
    pos   <- strsplit(sub("_meth\\.phen$", "", basename(file)), "_")[[1]]
    df    <- fread(meth_files[i], select = c("V1", "V3"))
    colnames(df) <- c("brnum", "meth")
    
    # Add in vmr pos
    df <- df %>%
      mutate(chr   = as.character(sub("chr_","", chr)),
             start = as.integer(pos[1]),
             end   = as.integer(pos[2]),
             feature_id = paste0("VMR", i))
    meth_list[[i]] <- df
  }
  # Bind meth matrix
  meth_df <- rbindlist(meth_list, use.names = TRUE, fill = TRUE)
  
  return(meth_df)
}

save_plot <- function(p, fn, w, h, dpi){
  for(ext in c(".pdf", ".png")){
    ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h, dpi=dpi)
  }
}

get_pca <- function(df, tissue, output_path) {
  filtered <- df %>%
    dplyr::select(brnum, meth, feature_id, h2_category)
  
  meth <- filtered %>%
    pivot_wider(names_from = brnum, values_from = meth) %>%
    column_to_rownames("feature_id")
  
  h2 <- filtered %>%
    distinct(feature_id, h2_category) %>%
    column_to_rownames("feature_id")
  
  meth_matrix <- as.matrix(meth %>% select(where(is.numeric)))
  
  pc <- prcomp(meth_matrix, scale = TRUE, center = TRUE)
  print(summary(pc))
  
  pc_df <- as.data.frame(pc$x)
  pc_df$h2_category <- h2[rownames(pc_df), "h2_category"]
  pc_df <- tibble::rownames_to_column(pc_df, "feature_id")
  
  # write file
  fwrite(pc_df, sep = "\t", file = file.path(output_path, 
                                             paste0(tissue, "_meth_pc.tsv")))

  return(pc_df)
}

plot_pca <- function(pc_df, h2_cols, output_path) {
  # format df for plotting 
  pc_df <- pc_df %>%
    select("PC1", "PC2", "h2_category")
  
  # plot pca
  p <- pc_df %>%
    ggplot(aes(x = PC1, y = PC2, 
               color = `h2_category`, 
               fill = `h2_category`)) +
    geom_point(size = 2, alpha = 0.4) +
    labs(x = "PC1",
         y = "PC2",
         color = "Heritability Category") +
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

tissues <- c("caudate", "dlpfc", "hippocampus")
h2_cols <- c(
  "Heritable" = "#497C8A",
  "Non-heritable" = "#8CA77B",
  "Low prediction" = "#E3A27F"
)

for (tissue in tissues) {
  # Get VMRs
  vmr <- get_vmrs(tissue)
  
  # Get methylation values
  meth_file_path <- here("heritability", tissue, "_m/vmr")
  meth_files     <- list.files(path = meth_file_path, pattern = "_meth\\.phen$", 
                               recursive = TRUE, full.names = TRUE)
  meth_df        <- merge_meth(meth_files)
  
  # Merge with h2 groups
  merged_df <- meth_df |>
    inner_join(vmr, by = c("chr" = "chrom", "start", "end"))
  
  # PCA
  pc_df  <- get_pca(merged_df, tissue, output_path)
  p      <- plot_pca(pc_df, h2_cols, output_path)
  pca_fn <- file.path(output_path, paste0(tissue, "_pca_meth"))
  save_plot(p, pca_fn, w = 10, h = 6, dpi = 300)
}

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()