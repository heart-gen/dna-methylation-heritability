## Generate circlized plot
suppressPackageStartupMessages({
    library(dplyr)
    library(here)
    library(circlize)
    library(data.table)
    library(ComplexHeatmap)
    library(RColorBrewer)
})

extract_bed <- function(num_indiv){
    target_file <- here("inputs", "simulated-data", "_m", 
                        paste0("sim_", num_indiv, "_indiv"), "snp_phenotype_mapping.tsv")
    target <- data.table::fread(target_file) %>%
      dplyr::select(chrom, start, end, target_heritability) %>%
      mutate(chrom = paste0("chr", chrom), 
             sample_size = paste0(num_indiv)) %>%
      rename(chr = chrom)
    
    target <- target %>%
      mutate(N = "Target", method = "Target",
             h2_category = case_when(
               target_heritability < 0.1 ~ "Non-heritable",
               target_heritability >= 0.1 ~ "Heritable"
             ))

    return(target)
}

plot_circos_7samples <- function(samples_list){
    circos.clear() # clear plot if there is any
    circos.par("start.degree" = 90, "gap.after" = c(rep(1, 23), 10),
              "track.height" = 0.06, "track.margin" = c(0.005, 0.005), 
              "cell.padding" = c(0, 0, 0, 0)) # rotate 90 degrees
                                        # initialize with ideogram
                                        # use hg38, default is hg19
    circos.initializeWithIdeogram(species="hg38")
    herit_colors <- c("Heritable" = "#0072B2", "Non-heritable" = "#1B9E77")
    
    bg_colors <- brewer.pal(n = length(samples_list), "Set3")
                                        # Plot each sample size as a track
    for (i in seq_along(samples_list)) {
      target <- samples_list[[i]]
      circos.genomicTrack(target, bg.col = add_transparency(bg_colors[i], transparency=0.8),
                          panel.fun = function(region, value, ...) {
        circos.genomicPoints(region, value, col = herit_colors[value$h2_category],
                             pch = 16, cex = 0.4)
      }, bg.border = bg_colors[i])
    }
    
                                        # Add stacked legends
    lgd_points <- Legend(at = names(herit_colors),
                  type = "points", pch = 16,
                  legend_gp = gpar(col = herit_colors, fontsize = 20),
                  title = "Heritability", background="#FFFFFF")
    
    lgd_tracks <- Legend(at = names(samples_list),
                  type = "box", pch = 16,
                  legend_gp = gpar(fill = add_transparency(bg_colors, transparency = 0.8),
                                   col = bg_colors, fontsize = 20),
                  title = "Sample size", background="#FFFFFF")
    
    lgd <- packLegend(lgd_points, lgd_tracks)
    
    draw(lgd, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))
}

output_path <- here("simulation-analysis/comparison/circos_plot/_m/")

if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

####### MAIN
main <- function(){
    sample_sizes <- c(100, 150, 200, 250, 500, 1000, 5000)
    
    samples_list <- lapply(sample_sizes, extract_bed)
    names(samples_list) <- paste0(sample_sizes, "_indiv")
    
                                        # plot
    pdf(file = file.path(output_path, "simulated_circos_plot.pdf"),
        width = 10, height = 10)
    plot_circos_7samples(samples_list)
    dev.off()
}

main()

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()