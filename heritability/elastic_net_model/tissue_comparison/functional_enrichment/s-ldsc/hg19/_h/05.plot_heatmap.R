library(ggplot2)
library(tidyverse)

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_annotation_enrichment <- function(file="results/ad/caudate/heritable_hg19/ad_caudate_heritable_hg19.results") {
  dt <- data.table::fread(file) %>%
    mutate(across(where(is.character), as.factor))
  return(dt)
}
memENRICH <- memoise::memoise(load_annotation_enrichment)

gen_data <- function() {
  err <- 1e-7
  dt <- memENRICH() %>%
    mutate(
      `-log10(P)` = -log10(Enrichment_p),
      `Enrichment Percentile` = Enrichment / (1 + Enrichment),
      p.sig = Enrichment_p < 0.05,
      `log2(Enrichment)` = log2(Enrichment + err),
      p.cat = cut(Enrichment_p, breaks = c(1, 0.05, 0.01, 0.005, 0),
                  labels = c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                  include.lowest = TRUE)
    )
  return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(label, w, h){
  df <- memDF() %>% filter(is.finite(`log2(Enrichment)`))
  
  y0 <- min(df$`log2(Enrichment)`, na.rm = TRUE) - 0.1
  y1 <- max(df$`log2(Enrichment)`, na.rm = TRUE) + 0.1
  
  tile_plot <- df %>%
    ggplot(aes(y = Category, x = Category, fill = `log2(Enrichment)`)) +
    geom_tile(color = "grey") +
    geom_text(aes(label = ifelse(p.sig,
                                 format(round(`-log10(P)`, 1), nsmall = 1), "")),
              color = "black", size = 5) +
    scale_fill_gradientn(colors = c("blue", "white", "red"),
                         values = scales::rescale(c(y0, 0, y1)),
                         limits = c(y0, y1),
                         name = "log2(Enrichment)") +
    labs(x = "Annotation", y = label) +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      axis.title = element_text(face = "bold", size = 28),
      strip.text = element_text(face = "bold", size = 22)
    )
  
  save_plot(tile_plot, paste0("tileplot_enrichment_", tolower(label)), w, h)
  print(tile_plot)
}
## Run script
plot_tile("Annotation Enrichment", 12, 12)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()