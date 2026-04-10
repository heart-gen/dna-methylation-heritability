library(ggplot2)
library(tidyverse)
library(tools)

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".png")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h, dpi=300)
    }
}

load_annotation_enrichment <- function(){
    return(data.table::fread("annotation_vmr_enrichment_analysis.txt"))
}
memENRICH <- memoise::memoise(load_annotation_enrichment)

## Annotation display labels (consistent order: genic first, then CpG, then distal)
ANNOT_LABELS <- c(
  "hg38_genes_promoters"  = "Promoter",
  "hg38_genes_1to5kb"     = "1\u20135 kb upstream",
  "hg38_genes_5UTRs"      = "5\u2032 UTR",
  "hg38_genes_3UTRs"      = "3\u2032 UTR",
  "hg38_genes_exons"      = "Exon",
  "hg38_genes_introns"    = "Intron",
  "hg38_enhancers_fantom" = "Enhancer",
  "hg38_cpg_islands"      = "CpG island",
  "hg38_cpg_shores"       = "CpG shore",
  "hg38_cpg_shelves"      = "CpG shelf",
  "hg38_cpg_inter"        = "Open sea",
  "hg38_genes_intergenic" = "Intergenic"
)

gen_data <- function(){
    err = 1e-7
    dt <- memENRICH() %>%
        filter(Annotation %in% names(ANNOT_LABELS)) %>%
        mutate(
            log2_or  = log2(OR + err),
            fdr_sig  = FDR < 0.05,
            sig_star = case_when(
                FDR < 0.001 ~ "***",
                FDR < 0.01  ~ "**",
                FDR < 0.05  ~ "*",
                TRUE        ~ ""
            ),
            Annotation = factor(ANNOT_LABELS[Annotation],
                                levels = rev(ANNOT_LABELS)),
            h2_Category = factor(h2_Category,
                                 levels = c("Heritable", "Non-heritable",
                                            "Low prediction")),
            Tissue = factor(Tissue,
                            levels = c("Caudate", "DLPFC", "Hippocampus"))
        )
    return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(label, w, h){
    df <- memDF() %>% filter(is.finite(log2_or))

    ## Symmetric colour limits so OR = 1 (log2 = 0) is always white
    lim <- max(abs(df$log2_or), na.rm = TRUE)
    lim <- ceiling(lim * 10) / 10          # round up to nearest 0.1

    tile_plot <- ggplot(df, aes(x = Tissue, y = Annotation, fill = log2_or)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = sig_star),
                  size = 4, vjust = 0.5, color = "black") +
        scale_fill_gradient2(
            low      = "#2166AC",
            mid      = "white",
            high     = "#B2182B",
            midpoint = 0,
            limits   = c(-lim, lim),
            name     = expression(log[2]~"(OR)"),
            breaks   = scales::pretty_breaks(n = 5)
        ) +
        facet_grid(. ~ h2_Category) +
        labs(x = NULL, y = NULL) +
        theme_classic(base_size = 11) +
        theme(
            axis.text.x        = element_text(size = 10),
            axis.text.y        = element_text(size = 10),
            strip.background   = element_blank(),
            strip.text         = element_text(face = "bold", size = 12),
            panel.spacing      = unit(0.8, "cm"),
            legend.position    = "right",
            legend.key.height  = unit(1.2, "cm"),
            legend.key.width   = unit(0.35, "cm"),
            legend.title       = element_text(size = 10),
            legend.text        = element_text(size = 9),
            plot.margin        = margin(6, 8, 6, 6)
        )

    save_plot(tile_plot, paste0("tileplot_enrichment_",tolower(label)), w, h)
}

## Run script
plot_tile("annotation", 8, 5.5)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()