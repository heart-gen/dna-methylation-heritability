#### Correlate VMR length to heritability estimates ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
})

## Function
summarise_h2 <- function(vmr, pop, out_path) {
  summary_df <- vmr %>%
    group_by(tissue) %>%
    summarise(
      n = n(),
      mean_h2 = mean(h2_unscaled, na.rm = TRUE),
      sd_h2 = sd(h2_unscaled, na.rm = TRUE),
      median_h2 = median(h2_unscaled, na.rm = TRUE),
      .groups = "drop"
    )
    write.csv(summary_df, file = file.path(out_path, paste0("h2_summary_all_sites_", pop, ".csv")), 
            row.names = FALSE)
  print(summary_df)
}

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}

plot_density <- function(vmr, tissue) {
  # Define palette
  tissue_colors <- c(
    "caudate" = "#7372A6",
    "dlpfc" = "#B36F61",
    "hippocampus" = "#C5AC47"
  )
  
  counts <- vmr %>%
    group_by(tissue) %>%
    summarise(n = n(), .groups = "drop")

  counts$tissue_label <- ifelse(
  tolower(counts$tissue) == "dlpfc", 
  "DLPFC", 
  tools::toTitleCase(counts$tissue)
  )

  labels <- setNames(
    paste0(counts$tissue_label, "\n(n = ", counts$n, ")"),
    counts$tissue
  )

  legend_labels <- setNames(counts$tissue_label, counts$tissue)
  
  p_hist <- gghistogram(vmr, x = "h2_unscaled", 
                        add_density = TRUE, rug = TRUE, 
                        add = "mean",
                        color = "tissue", fill = "tissue",
                        ggtheme = theme_pubr(base_size = 20, border = TRUE),
                        xlab = "Estimated h²", ylab = "Count") +
    facet_wrap(~tissue, labeller = as_labeller(labels), scales = "free_x") +
    scale_color_manual(values = tissue_colors, labels = legend_labels) +
    scale_fill_manual(values = tissue_colors, labels = legend_labels) +
    labs(color = NULL, fill = NULL) +
    font("xy.title", face = "bold", size = 14) +
    geom_vline(xintercept = 0.19, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0.1, linetype = "dotted", 
               color = "#8CA77B", size  = 2) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  return(p_hist)
}

## Main
tissues <- c("caudate", "hippocampus", "dlpfc")
pops <- c("AA", "EA")

out_path <- here("heritability/elastic_net_model/all_individuals/tissue_comparison/h2_distribution/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# Get all VMRs per population
vmr_all_pop <- list()

for (pop in pops) {

  vmr_all_tissue <- list()

  for (tissue in tissues) {
    # Read in summary table
    enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                      paste0(tissue, "/_m/", tissue, "_summary_elastic-net_", 
                      pop, ".tsv"))
    vmr <- read.table(enet_file, sep = "\t", header = TRUE, 
                      colClasses = c(chrom = "character")) |> na.omit()
      
    # Store vmrs across all tissues
    vmr$tissue <- tissue
    vmr_all_tissue[[tissue]] <- vmr
  }

  vmr_all_pop[[pop]] <- bind_rows(vmr_all_tissue)

}

vmr_all_AA <- vmr_all_pop[["AA"]]
vmr_all_EA <- vmr_all_pop[["EA"]]

# Get matched VMRs across populations
vmr_matched    <- list()

for (tissue in tissues) {
    # Read in summary table
    enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                      paste0(tissue, "/_m/", tissue, "_summary_elastic-net_matched_r2_0.3.tsv"))
    vmr <- read.table(enet_file, sep = "\t", header = TRUE, 
                      colClasses = c(chrom = "character")) |> na.omit()
      
    # Store vmrs across all tissues
    vmr$tissue <- tissue
    vmr_matched[[tissue]] <- vmr
}

vmr_matched <- bind_rows(vmr_matched)

vmr_matched_AA <- vmr_matched %>% 
    select(chrom, start, end, tissue, h2_unscaled = h2_unscaled_AA, h2_category)

vmr_matched_EA <- vmr_matched %>% 
    select(chrom, start, end, tissue, h2_unscaled = h2_unscaled_EA, h2_category)

# Define vmr groups to loop through
vmr_groups <- list(
  AA_all = vmr_all_AA,
  EA_all = vmr_all_EA,
  AA_matched = vmr_matched_AA,
  EA_matched = vmr_matched_EA
)

for (group in names(vmr_groups)){

  vmr <- vmr_groups[[group]]

  # Summarize h2 of vmrs
  summarise_h2(vmr, group, out_path)

  # Save plot
  p_hist <- plot_density(vmr, tissue)
  fn_hist <- file.path(out_path, paste0("all_sites_vmr_h2_distribution_", group))
  save_plot(p_hist, fn_hist, 10, 5)
}

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
