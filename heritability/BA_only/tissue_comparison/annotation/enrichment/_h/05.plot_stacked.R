#### Plot top enrichment across brain regions ####

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(data.table)
  library(tidyverse)
  library(ggplot2)
})

## Function 
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

plot_enrichment <- function(df, tissue_cols, h2_cat){
  top_enrichment <- df |> 
    filter(h2_Category == h2_cat, FDR < 0.05) |>
    group_by(Tissue) |>
    arrange(`-log10(FDR)`) |> slice_head(n = 5) |> ungroup()

  p = ggplot(top_enrichment, aes(x=`-log10(FDR)`, y = Annotation,
                                 color = Tissue, size=`log2(OR)`)) +
    geom_point(shape=18, alpha=0.8) +
    labs(y='', x='-Log10 (Adjusted P-value)', size="Enrichment") +
    scale_colour_manual(name="Brain Region", values=tissue_cols,
                        labels=c("Caudate","DLPFC","Hippocampus")) +
    scale_size_continuous(range = c(2, 10)) +
    guides(
      colour = guide_legend(override.aes = list(size = 6))
    ) +
    theme_bw(base_size=15) +
    theme(axis.title=element_text(face='bold'),
          strip.text=element_text(face='bold'))
  
  return(p)
}

## Main

# Load annotation results
annot <- memoise::memoise(load_annotation_enrichment)
memDF <- memoise::memoise(gen_data)
df <- memDF() %>% filter(is.finite(`log2(OR)`))

# Create output dir
output_path <- here("heritability", "elastic_net_model", "BA_only",
                    "tissue_comparison", "annotation", "enrichment", "_m",
                    "stacked_plot")

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

h2_categories <- c("Heritable", "Non-heritable", "Low prediction")
tissue_cols <- c(
  "Caudate" = "#B36F61",
  "DLPFC" = "#7372A6",
  "Hippocampus" = "#E3C962"
)

for (h2_cat in h2_categories) {
  # Plot enrichment across brain regions
  p <- plot_enrichment(df, tissue_cols, h2_cat)
  stacked_fn <- file.path(output_path, paste0(tolower(h2_cat), "_stacked"))
  save_plot(p, stacked_fn, w = 10, h = 6, dpi = 300)
}

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()