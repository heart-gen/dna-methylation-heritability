##### Generates plots for GO enrichment analysis. #####
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
    library(here)
})

# Function
save_plot <- function(p, fn, w, h){
    for(ext in c('.png', '.pdf')){
        ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

get_top_GO <- function(tissue, gene_set){
  err <- 1e-15
  fn  <- here(file.path("heritability/elastic_net_model/tissue_comparison/functional_enrichment/_m", 
                        paste0(tolower(tissue), "_heritable_h2_r2_0.75_", gene_set, ".csv")))
  return(data.table::fread(fn) |>
           filter(stringr::str_detect(id, "^GO")) |>
           arrange(p_value) |> head(10) |>
           mutate(`Log10`=-log10(p_adjust+err), Tissue=tissue))
}

get_top_KEGG <- function(tissue){
  err <- 1e-15
  fn  <- here(file.path("heritability/elastic_net_model/tissue_comparison/functional_enrichment/_m", 
                        paste0(tolower(tissue), "_heritable_h2_r2_0.75_KEGG.csv")))
  
  # Load KEGG pathway descriptions
  kegg_map <- read.delim("https://rest.kegg.jp/list/pathway/hsa", 
                         header = FALSE, sep = "\t")
  colnames(kegg_map) <- c("id", "description")
  kegg_map$id <- sub("path:", "", kegg_map$id)
  
  # Read KEGG enrichment results
  dt <- data.table::fread(fn)
  
  # Ensure 'id' is character type to match KEGG map
  dt$id <- as.character(dt$id)
  
  # Get top KEGG 
  dt <- dt |>
    arrange(p_value) |> head(10) |>
    mutate(Log10 = -log10(p_adjust + err), Tissue = tissue)
  
  # Merge in descriptions
  dt <- left_join(dt, kegg_map, by = "id")
  
  return(dt)
}

generate_dataframe <- function(gene_set){
    df_list <- list()
    tissues <- c("Caudate", "DLPFC", "Hippocampus")
    
    for(jj in seq_along(tissues)){
      if (gene_set %in% c("GO_BP", "GO_MF")){
        df_list[[jj]] <- get_top_GO(tissues[jj], gene_set)
      } else if (gene_set == "KEGG"){
        df_list[[jj]] <- get_top_KEGG(tissues[jj])
      } else {
        stop("Invalid gene set")
      }
    }
    return( bind_rows(df_list) )
}

plot_enrichment <- function(gene_set){
    dt <- generate_dataframe(gene_set)
    cbPalette <- brewer.pal(4, "Set1")
    gg1 = ggplot(dt, aes(x=Log10, y=description, color=Tissue,
                         size=fold_enrichment)) +
        geom_point(shape=18, alpha=0.8) +
        labs(y='', x='-Log10 (Adjusted P-value)', size="Enrichment") +
        scale_colour_manual(name="Brain Region", values=cbPalette,
                            labels=c("Caudate","DLPFC","Hippocampus")) +
        scale_size_continuous(range = c(2, 10)) +
        guides(
          colour = guide_legend(override.aes = list(size = 6))
        ) +
        theme_bw(base_size=15) +
        theme(axis.title=element_text(face='bold'),
              strip.text=element_text(face='bold'))
    return(gg1)
}

# Main
out_path = here("heritability/elastic_net_model/tissue_comparison/functional_enrichment/_m/plots")
if (!dir.exists(out_path)) {
        dir.create(out_path, recursive = TRUE)
}

for (gene_set in c("GO_BP", "GO_MF", "KEGG")){
  gg = plot_enrichment(gene_set)
  fn = file.path(out_path, paste0("heritable_VMRs_r2_0.75.", gene_set, ".stacked"))
  save_plot(gg, fn, 14, 6)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
