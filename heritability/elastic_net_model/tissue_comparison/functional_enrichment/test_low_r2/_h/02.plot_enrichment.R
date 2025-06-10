##### Generates plots for GO enrichment analysis. #####
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(RColorBrewer)
})

# Function
save_plot <- function(p, fn, w, h){
    for(ext in c('.png', '.pdf')){
        ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

get_top_GO <- function(tissue){
    err <- 1e-15
    fn  <- file.path("/projects/p32505/users/alexis/working/_m", 
                paste0(tolower(tissue), "_heritable_h2_r2_0.75_GO_BP.csv"))
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
        guides(
          colour = guide_legend(override.aes = list(size = 6))
        ) +
        theme_bw(base_size=15) +
        theme(axis.title=element_text(face='bold'),
              strip.text=element_text(face='bold'))
    return(gg1)
}

# Main
gg = plot_GO()
out_path = file.path("/projects/p32505/users/alexis/working/_m/plots")
if (!dir.exists(out_path)) {
        dir.create(out_path, recursive = TRUE)
}
fn = file.path(out_path, "heritable_VMRs_r2_0.75.GO_BP.stacked")
save_plot(gg, fn, 14, 6)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
