suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(ggpubr)
  library(txdbmaker)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(data.table)
})


get_dist <- function(annotation) {

  vmr_gr <- GRanges(
    seqnames = annotation$seqnames,
    ranges = IRanges(start = annotation$start, end = annotation$end)
  )
  
  genes_gr <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  nearest_gene <- distanceToNearest(vmr_gr, genes_gr)
  
  annotation <- annotation %>%
    mutate(distance_to_nearest_gene = mcols(nearest_gene)$distance)
  
  return(annotation)
}

summarise_dist <- function(annotation, tissue, out_path) {
  summary_df <- annotation %>%
    group_by(h2_category, annot.type) %>%
    summarise(
      n = n(),
      mean_dist = mean(distance_to_nearest_gene, na.rm = TRUE),
      sd_dist = sd(distance_to_nearest_gene, na.rm = TRUE),
      median_dist = median(distance_to_nearest_gene, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_df)
  write.csv(summary_df, file = file.path(out_path, paste0("distance_summary_", tissue, ".csv")), 
            row.names = FALSE)
}

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}

plot_density <- function(annotation, tissue) {
  category_colors <- c(
    "caudate"     = "#B36F61", 
    "dlpfc"       = "#7372A6", 
    "hippocampus" = "#E3C962"  
  )
  
  tissue_title <- ifelse(tolower(tissue) == "dlpfc", "DLPFC", tools::toTitleCase(tissue))
  
  p_hist <- gghistogram(annotation, x = "distance_to_nearest_gene", 
                        add_density = TRUE, rug = TRUE, 
                        add = "median",
                        ggtheme = theme_pubr(base_size = 15, border = TRUE),
                        xlab = "Distance to nearest gene (bp)", ylab = "Count") +
    scale_color_manual(values = category_colors) +
    scale_fill_manual(values = category_colors) +
    ggtitle(paste("Distance to nearest gene:", tissue_title)) +
    labs(color = NULL, fill = NULL) +
    font("xy.title", face = "bold", size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p_hist)
}

# Main
tissues <- c("caudate", "hippocampus", "dlpfc")
out_path <- here("heritability/elastic_net_model/BA_only/tissue_comparison/annotation/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

hist_plots <- list()

for (tissue in tissues) {
  # Read in summary table
  annotation_file <- here("heritability/elastic_net_model/BA_only/tissue_comparison/annotation/_m", 
                    paste0(tissue, "_vmr_annotations_hg38.tsv"))
  annotation <- fread(annotation_file, sep = "\t")
  
  # Add distance to annot 
  annotation <- get_dist(annotation)
  
  # Summarize
  summarise_dist(annotation, tissue, out_path)
  
  # Plotting
  p_hist <- plot_density(annotation, tissue)
  
  # Store plots
  hist_plots[[tissue]] <- p_hist
}

# Save plots
combined_hist <- ggarrange(plotlist = hist_plots, ncol = 3, nrow = 1)
fn_hist <- file.path(out_path, "gene_dist_distribution")
save_plot(combined_hist, fn_hist, 14, 6)
