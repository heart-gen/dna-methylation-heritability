#### Correlate enrichment across brain regions ####

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
  library(data.table)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
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

save_plot <- function(p, fn, w, h){
  for(ext in c(".pdf", ".png")){
    ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
  }
}

plot_scatter <- function(df, h2_cat, tissue1, tissue2, output_path){
  pivot <- df %>%
    filter(h2_Category == h2_cat) %>%
    select(Tissue, Annotation, `log2(OR)`) %>%
    pivot_wider(names_from = Tissue, values_from = `log2(OR)`)
  
  xlab = paste0("log2(OR)\n", tissue1)
  ylab = paste0("log2(OR)\n", tissue2)
  
  fn = file.path(output_path, paste(tolower(h2_cat), "effectsize_scatter", 
             tolower(tissue1), tolower(tissue2), sep="_"))
  
  pp <- ggscatter(pivot, x = tissue1, y = tissue2, add = "reg.line", 
                  size = 1, xlab = xlab, ylab = ylab, panel.labs.font=list(face="bold"),
                  add.params=list(color = "blue", fill = "lightgray"),
                  conf.int = TRUE, cor.coef = TRUE, cor.coef.size = 3,
                  cor.method = "spearman", cor.coeff.args=list(label.sep="\n"), ncol = 4) +
    theme_pubr(base_size=18)
  save_plot(pp, fn, 6, 6)
}

# Main
output_path <- here("heritability", "elastic_net_model", "tissue_comparison",
                    "annotation", "enrichment", "_m", "scatter_plot")
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

annot <- memoise::memoise(load_annotation_enrichment)
memDF <- memoise::memoise(gen_data)
df <- memDF() %>% filter(is.finite(`log2(OR)`))

for (h2_cat in c("Heritable", "Non-heritable", "Low prediction")){
  for (tissue1 in c("Caudate", "DLPFC", "Hippocampus")){
    for (tissue2 in c("DLPFC", "Hippocampus")){
      if (tissue1 != tissue2){
        plot_scatter(df, h2_cat, tissue1, tissue2, output_path)
      }
    }
  }
}

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()