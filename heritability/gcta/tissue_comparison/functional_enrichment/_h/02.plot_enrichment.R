##### Generates plots for GO enrichment analysis. #####
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(here)
    library(RColorBrewer)
})

save_plot <- function(p, fn, w, h){
    for(ext in c('.png', '.pdf')){
        ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

get_top_GO <- function(tissue){
    err <- 1e-15
    fn  <- here("heritability/gcta/tissue_comparison/functional_enrichment/_m", 
                paste0(tolower(tissue), "_non_heritable_h2_pval_GO.csv"))
    return(data.table::fread(fn) |>
           filter(stringr::str_detect(id, "^GO")) |>
           arrange(p_value) |> head(10) |>
           mutate(`Log10`=-log10(p_adjust+err), Tissue=tissue))
}

generate_dataframe <- function(){
    df_list <- list()
    tissues <- c("Caudate", "DLPFC", "Hippocampus")
    for(jj in seq_along(tissues)){
        df_list[[jj]] <- get_top_GO(tissues[jj])
    }
    return( bind_rows(df_list) )
}

plot_GO <- function(){
    dt <- generate_dataframe()
    cbPalette <- brewer.pal(4, "Set1")
    gg1 = ggplot(dt, aes(x=Log10, y=description, color=Tissue,
                         size=fold_enrichment)) +
        geom_point(shape=18, alpha=0.8) +
        labs(y='', x='-Log10 (Adjusted P-value)', size="Enrichment") +
        scale_colour_manual(name="Brain Region", values=cbPalette,
                            labels=c("Caudate","DLPFC","Hippocampus")) +
        scale_size_continuous(range = c(2, 10)) +
        theme_bw(base_size=15) +
        theme(axis.title=element_text(face='bold'),
              strip.text=element_text(face='bold'))
    return(gg1)
}

#### MAIN
gg = plot_GO()
fn = here("heritability/gcta/tissue_comparison/functional_enrichment/_m/non_heritable_VMRs.GO.stacked")
save_plot(gg, fn, 14, 6)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
