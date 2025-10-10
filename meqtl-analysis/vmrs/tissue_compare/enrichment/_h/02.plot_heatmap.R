library(ggplot2)
library(tidyverse)

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_meqtl_enrichment <- function(){
    return(data.table::fread("meqtl_vmr_enrichment_analysis.txt"))
}
memENRICH <- memoise::memoise(load_meqtl_enrichment)

gen_data <- function(){
    err = 0.0000001
    dt <- memENRICH() %>% mutate(across(where(is.character), as.factor)) %>%
        mutate(h2_Category=fct_relevel(h2_Category, rev), `-log10(FDR)`= -log10(FDR),
               `OR Percentile`= OR / (1+OR), p.fdr.sig=FDR < 0.05,
               `log2(OR)` = log2(OR+err),
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE)) %>%
        mutate(across(Direction, factor, levels=c("All", "Down", "Up")))
    #levels(dt$Direction) <- c("All", "Heritable Bias", "Non-heritable Bias")
    return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(label, w, h){
  df <- memDF() %>% filter(Direction == "All", is.finite(`log2(OR)`))
  
  y0 <- min(df$`log2(OR)`, na.rm = TRUE) - 0.1
  y1 <- max(df$`log2(OR)`, na.rm = TRUE) + 0.1
  
  tile_plot <- df %>%
    ggplot(aes(y = h2_Category, x = Tissue, fill = `log2(OR)`)) +
    geom_tile(color = "grey") +
    geom_text(aes(label = ifelse(p.fdr.sig,
                                 format(round(`-log10(FDR)`,1), nsmall=1), "")),
              color = "black", size = 5) +
    scale_fill_gradientn(colors = c("blue", "white", "red"),
                         values = scales::rescale(c(y0, 0, y1)),
                         limits = c(y0, y1),
                         name = "log2(OR)") +
    #facet_grid(. ~ Direction) +
    labs(x = "Tissue", y = label) +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      axis.title = element_text(face = "bold", size = 28),
      strip.text = element_text(face = "bold", size = 22)
    )
  
  save_plot(tile_plot, paste0("tileplot_enrichment_",tolower(label)), w, h)
  print(tile_plot)
}
## Run script
plot_tile("Heritability Category", 10, 10)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()