## Generate circlized plot
suppressPackageStartupMessages({
    library(dplyr)
    library(here)
    library(circlize)
    library(ComplexHeatmap)
})

extract_bed <- function(tissue){
    enet_file <- here("heritability/elastic_net_model/", 
                      paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
    vmr <- data.table::fread(enet_file) %>%
      dplyr::select(chrom, start, end, h2_unscaled, r_squared_cv) %>%
      mutate(chrom = paste0("chr", chrom)) %>%
      rename(chr = chrom)
    
    vmr <- vmr %>%
      mutate(h2_category = case_when(
        r_squared_cv <= 0.75 ~ "Low prediction",
        h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
        h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
      ))
    
    vmr_groups <- split(vmr, vmr$h2_category)
    
    return(vmr_groups)
}

plot_circos_3tissue <- function(caudate, dlpfc, hippo){
    lgd_points = Legend(at=c("Heritable", "Non-heritable", "Low prediction"),
                        type="points", legend_gp=gpar(col = c("red", "blue", "green")),
                        title_position="topleft", title="Heritability",
                        background="#FFFFFF")
    circos.clear() # clear plot if there is any
    circos.par("start.degree" = 0, "cell.padding" = c(0, 0, 0, 0),
               "track.height" = 0.15) # rotate 90 degrees
    # initialize with ideogram
    # use hg38, default is hg19
    circos.initializeWithIdeogram(species="hg38")
    circos.genomicTrack(caudate, bg.border="#E64B35FF",
                        bg.col=add_transparency("#E64B35FF", transparency=0.7),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.6, col = c("red", "blue", "green")[i], ...)
    })
    circos.genomicTrack(dlpfc, bg.border="#00A087FF",
                        bg.col=add_transparency("#00A087FF", transparency=0.7),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.6, col = c("red", "blue", "green")[i], ...)
    })
    circos.genomicTrack(hippo, bg.border="#3C5488FF",
                        bg.col=add_transparency("#3C5488FF", transparency=0.7),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.6, col = c("red", "blue", "green")[i], ...)
    })
    draw(lgd_points, x=unit(5, "mm"), y=unit(5, "mm"), just=c("left", "bottom"))
}

####### MAIN
main <- function(){
    caudate <- extract_bed("caudate")
    dlpfc   <- extract_bed("dlpfc")
    hippo   <- extract_bed("hippocampus")
                                        # plot
    pdf(file = paste0("significant_circos_plot_3regions.pdf"),
        width = 10, height = 10)
    plot_circos_3tissue(caudate, dlpfc, hippo)
    dev.off()
}

main()

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
