## Generate circlized plot for VMRs
suppressPackageStartupMessages({
    library(dplyr)
    library(here)
    library(circlize)
    library(ComplexHeatmap)
})

extract_bed <- function(tissue){
    enet_file <- here("heritability/elastic_net_model/BA_only/", 
                      paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
    vmr <- data.table::fread(enet_file) %>%
      mutate(chrom = paste0("chr", chrom)) %>%
      rename(chr = chrom) %>%
      na.omit() %>%
      mutate(h2_category = case_when(
        r_squared_cv <= 0.3 ~ "Low prediction",
        h2_unscaled < 0.1 & r_squared_cv > 0.3 ~ "Non-heritable",
        h2_unscaled >= 0.1 & r_squared_cv > 0.3 ~ "Heritable"
      ))
    
    vmr$start <- vmr$start - 500000
    vmr$end   <- vmr$end + 500000
    vmr       <- vmr %>% select(chr, start, end, h2_unscaled, h2_category) # keep for plotting
    
    vmr_groups <- split(vmr, vmr$h2_category)
    
    return(vmr_groups)
}

plot_circos_3tissue <- function(caudate, dlpfc, hippo, col, ylim){
                                            # set custom gap for y axis labels
    chrs <- paste0("chr", c(1:22, "X", "Y"))
    gap <- rep(1, length(chrs))
    names(gap) <- chrs
    gap["chrY"] <- 9
  
    lgd_points = Legend(at=c("Caudate nucleus", "DLPFC", "Hippocampus"),
                        type="points", legend_gp=gpar(col = c("#7372A6", "#B36F61", "#C5AC47")),
                        title_position="topleft", title="Tissue",
                        background="#FFFFFF")
    lgd_threshold <- Legend(
      labels = expression("h"^2*" = 0.1 threshold"),
      type = "lines",
      legend_gp = gpar(lty = 2, col = "#6E6E6E"),
      background = "#FFFFFF"
    )
    circos.clear() # clear plot if there is any
    circos.par("start.degree" = 90, "cell.padding" = c(0, 0, 0, 0),
               "track.height" = 0.15, gap.degree = gap) # rotate 90 degrees
    # initialize with ideogram
    # use hg38, default is hg19
    circos.initializeWithIdeogram(species="hg38")
    circos.genomicTrack(caudate, ylim = ylim, bg.border="#7372A6",
                        bg.col=add_transparency("#7372A6", transparency=0.8),
                        panel.fun = function(region, value, ...) {
                          circos.genomicLines(region, value,
                                              type = "segment", col = col, lwd = 2, ...)
                          circos.lines(CELL_META$cell.xlim, c(0.1, 0.1), lty = 2, col = "black")
                  			  circos.yaxis(labels.cex = 0.6, side = "left", sector.index = "chr1")
    })
    circos.genomicTrack(dlpfc, ylim = ylim, bg.border="#B36F61",
                        bg.col=add_transparency("#B36F61", transparency=0.8),
                        panel.fun = function(region, value, ...) {
                          circos.genomicLines(region, value,
                                              type = "segment", col = col, lwd = 2, ...)
                          circos.lines(CELL_META$cell.xlim, c(0.1, 0.1), lty = 2, col = "black")
    })
    circos.genomicTrack(hippo, ylim = ylim, bg.border="#C5AC47",
                        bg.col=add_transparency("#C5AC47", transparency=0.8),
                        panel.fun = function(region, value, ...) {
                          circos.genomicLines(region, value,
                                              type = "segment", col = col, lwd = 2, ...)
                          circos.lines(CELL_META$cell.xlim, c(0.1, 0.1), lty = 2, col = "black")
    })
    lgd_combined <- packLegend(lgd_points, lgd_threshold, direction = "vertical")
    draw(lgd_combined, x=unit(5, "mm"), y=unit(5, "mm"), just=c("left", "bottom"))
}

                                        # create output dir
output_path <- here("heritability/elastic_net_model/BA_only/tissue_comparison/circos_plot/_m/")

if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

####### MAIN
main <- function(){
                                        # extract filtered dataframes
    caudate <- extract_bed("caudate") 
    dlpfc   <- extract_bed("dlpfc")
    hippo   <- extract_bed("hippocampus")
    
                                        # define h2 categories and params
    categories <- c("Heritable", "Non-heritable", "Low prediction")
    category_colors <- list(
      "Heritable"      = "#6E6E6E",
      "Non-heritable"  = "#6E6E6E",
      "Low prediction" = "#6E6E6E"
    )
    category_ylims <- list(
      "Heritable"      = c(0, 1),
      "Non-heritable"  = c(0, 0.15),
      "Low prediction" = c(0, 0.6)
    )
    
                                        # plot circos for each h2 category
    for (category in categories) {
      caudate_cat <- caudate[[category]]
      dlpfc_cat   <- dlpfc[[category]]
      hippo_cat   <- hippo[[category]]
      
      pdf(file = file.path(output_path, 
                           paste0("free_y_circos_plot_3regions_", gsub(" ", "_", tolower(category)), ".pdf")),
          width = 5, height = 5)
      plot_circos_3tissue(caudate_cat, dlpfc_cat, hippo_cat, 
                          col  = category_colors[[category]],
                          ylim = category_ylims[[category]])
      dev.off()
    }
}

main()

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
