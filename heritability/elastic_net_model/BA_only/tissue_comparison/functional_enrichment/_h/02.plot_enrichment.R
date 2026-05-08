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

get_top_GO <- function(tissue, h2_cat, gene_set, top_n = 10){
  err <- 1e-15
  fn  <- here(file.path("heritability/elastic_net_model/BA_only/tissue_comparison/functional_enrichment/_m", 
                        gene_set, paste0(tolower(tissue), "_", h2_cat, ".csv")))
  return(data.table::fread(fn) |>
           filter(stringr::str_detect(id, "^GO")) |>
           filter(p_adjust < 0.05) |>
           arrange(p_adjust) |> head(top_n) |>
           mutate(`Log10`=-log10(p_adjust+err), Tissue=tissue))
}

get_top_KEGG <- function(tissue, h2_cat, top_n = 10){
  err <- 1e-15
  fn  <- here(file.path("heritability/elastic_net_model/BA_only/tissue_comparison/functional_enrichment/_m/KEGG", 
                        paste0(tolower(tissue), "_", h2_cat, ".csv")))
  
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
    filter(p_adjust < 0.05) |>
    arrange(p_adjust) |> head(top_n) |>
    mutate(Log10 = -log10(p_adjust + err), Tissue = tissue)
  
  # Merge in descriptions
  dt <- left_join(dt, kegg_map, by = "id")
  
  return(dt)
}

generate_dataframe <- function(gene_set, h2_cat){
    df_list <- list()
    tissues <- c("Caudate", "DLPFC", "Hippocampus")
    
    for(jj in seq_along(tissues)){
      if (gene_set %in% c("GO_BP", "GO_MF")){
        df_list[[jj]] <- get_top_GO(tissues[jj], h2_cat, gene_set)
      } else if (gene_set == "KEGG"){
        df_list[[jj]] <- get_top_KEGG(tissues[jj], h2_cat)
      } else if (gene_set == "GO_BP_KEGG"){
        df_list[[jj]] <- bind_rows(
          get_top_KEGG(tissues[jj], h2_cat, top_n = 3),
          get_top_GO(tissues[jj], h2_cat, "GO_BP", top_n = 3))
      } else{
        stop("Invalid gene set")
      }
    }
    return( bind_rows(df_list) )
}

plot_enrichment <- function(gene_set, h2_cat){
    dt <- generate_dataframe(gene_set, h2_cat)
    tissue_cols <- c(
      "Caudate" = "#B36F61",
      "DLPFC" = "#7372A6",
      "Hippocampus" = "#E3C962"
    )
    gg1 = ggplot(dt, aes(x=Log10, y=description, color=Tissue,
                         size=fold_enrichment)) +
        geom_point(shape=18, alpha=0.8) +
        labs(y='', x='-Log10 (Adjusted P-value)', size="Enrichment") +
        scale_colour_manual(name="Brain Region", values=tissue_cols) +
        scale_size_continuous(range = c(2, 10)) +
        guides(
          colour = guide_legend(override.aes = list(size = 6))
        ) +
        theme_bw(base_size=20) +
        theme(axis.title=element_text(face='bold'),
              strip.text=element_text(face='bold'))
    return(gg1)
}

# Main
for (h2_cat in c("heritable", "non_heritable", "low_prediction", "all")){
  for (gene_set in c("GO_BP", "GO_MF", "KEGG", "GO_BP_KEGG")){
    # create subdirs for plots
    out_path = here("heritability/elastic_net_model/BA_only/tissue_comparison/functional_enrichment/_m/", gene_set, "plots")
    if (!dir.exists(out_path)) {
      dir.create(out_path, recursive = TRUE)
    }
    
    # plot enrichment
    gg = plot_enrichment(gene_set, h2_cat)
    fn = file.path(out_path, paste0(h2_cat, "_VMRs", ".stacked"))
    save_plot(gg, fn, 14, 8)
  }
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
