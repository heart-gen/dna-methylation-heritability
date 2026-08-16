# Visualize h2 and r2 thresholds for classification

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

## Function
filter_sites <- function(enet) {
  vmr <- na.omit(enet)
  vmr <- vmr %>%
    mutate(h2_category = case_when(
      r_squared_cv <= 0.3 ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.3 ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.3 ~ "Heritable"
    ),
    h2_category = factor(h2_category, levels = c("Heritable", 
                                                 "Non-heritable", 
                                                 "Low prediction"))
    )
  return(vmr)
}

plot_scatter <- function(vmr) {
  p <- ggplot(vmr, aes(x = h2_unscaled, y = r_squared_cv, color = h2_category)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  geom_hline(yintercept = 0.3, linetype = "dashed") +
  scale_color_manual(
    values = c("Heritable" = "#497C8A",
               "Non-heritable" = "#8CA77B",
               "Low prediction" = "#E3A27F")
  ) +
  labs(
    x = "h2",
    y = "r2",
    color = "Category"
  ) +
  facet_wrap(~tissue) +
  theme_minimal(base_size = 20) + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.title = element_text(hjust = 0.5),
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

  return(p)
}

save_plot <- function(p, fn, w, h, dpi){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn,ext), plot=p, width=w, height=h, dpi=dpi)
  }
}

## Main 
# Create output directory
out_path <- here("heritability/elastic_net_model/all_individuals/tissue_comparison/h2_distribution/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# Read in and filter summary tables
tissues <- c("caudate", "hippocampus", "dlpfc")
pops <- c("AA", "EA")

# Get all VMRs per population
vmr_all_pop <- list()

for (pop in pops) {

  vmr_all_tissue <- list()

  for (tissue in tissues) {
    # Read in summary table
    enet_file <- here("heritability/elastic_net_model/all_individuals/", 
                      paste0(tissue, "/_m/", tissue, "_summary_elastic-net_", 
                      pop, ".tsv"))
    enet <- read.table(enet_file, sep = "\t", header = TRUE, 
                      colClasses = c(chrom = "character"))

    vmr <- filter_sites(enet)
      
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
    select(chrom, start, end, tissue, r_squared_cv = r_squared_cv_EA,     
           h2_unscaled = h2_unscaled_AA, h2_category)

vmr_matched_EA <- vmr_matched %>% 
    select(chrom, start, end, tissue, r_squared_cv = r_squared_cv_EA,
           h2_unscaled = h2_unscaled_EA, h2_category)

# Define vmr groups to loop through
vmr_groups <- list(
  AA_all = vmr_all_AA,
  EA_all = vmr_all_EA,
  AA_matched = vmr_matched_AA,
  EA_matched = vmr_matched_EA
)

for (group in names(vmr_groups)){

  vmr <- vmr_groups[[group]]

  # Plot
  p <- plot_scatter(vmr)

  # Save plot
  plot_file <- file.path(out_path, paste0("h2_r2_scatter_", group))
  save_plot(p, plot_file, w = 8, h = 4, dpi = 300)

}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
