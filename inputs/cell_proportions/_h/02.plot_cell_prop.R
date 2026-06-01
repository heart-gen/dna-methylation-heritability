#### Plot cell proportions from deconvolution ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
})

## Function
clean_pheno <- function(pheno_file_path, samples_to_include){
  pheno_df <- fread(pheno_file_path, header = TRUE) |>
    inner_join(samples_to_include, by = c("brnum", "region")) |>
    filter(agedeath > 17) |>
    mutate(region = recode(region,
                           "caudate" = "Caudate",
                           "dlpfc" = "DLPFC",
                           "hippocampus" = "Hippocampus"),
            snpPC1 = as.numeric(snpPC1),
            snpPC2 = as.numeric(snpPC2),
            donor_group = recode(race, "AA" = "BA", "EA" = "WA")) |>
    mutate_if(is.character, as.factor)
    
  return(pheno_df)
}

summarize_prop <- function(cell_prop, region, out_path){
    summary_df <- cell_prop %>%
    group_by(cell_type) %>%
    summarise(
      n = n(),
      mean = mean(proportion, na.rm = TRUE),
      median = median(proportion, na.rm = TRUE),
      sd = sd(proportion, na.rm = TRUE),
      min = min(proportion, na.rm = TRUE),
      max = max(proportion, na.rm = TRUE),
      .groups = "drop"
    )
  print(summary_df)
  write.csv(summary_df, 
            file = file.path(out_path, 
                             paste0("music-proportions-summary-", region, ".csv")), 
            row.names = FALSE)
}

plot_stack <- function(cell_prop){

    cell_type_colors <- c(
    "Astro"  = "#E64B35", 
    "Excit"  = "#4DBBD5", 
    "Immune" = "#91D1C2",
    "Inhib"  = "#00A087", 
    "Micro"  = "#6B6ECF", 
    "D1-SPN" = "#7E6148",
    "D2-SPN" = "#B09C85", 
    "Mural"  = "#DC91C0", 
    "Oligo"  = "#3C5488",
    "OPC"    = "#F39B7F"
    )
    celltype_order <- names(cell_type_colors)

    sample_order <- cell_prop |>
        select(brnum, donor_group) |>
        distinct() |>
        arrange(donor_group, brnum) |>
        pull(brnum)

    cell_prop <- cell_prop |>
        mutate(brnum     = factor(brnum, levels = sample_order),
               cell_type = factor(cell_type, levels = celltype_order))

    donor_group_colors <- c(
    "BA" = "#a52a2a66",
    "WA" = "#0000ff66"
    )

    donor_group_annot <- cell_prop |>
        select(brnum, donor_group) |>
        distinct()

    ## Main stacked barplot
    p_main <- ggplot(cell_prop, aes(x = brnum, y = proportion, fill = cell_type)) +
        geom_col(position = "stack", width = 0.75) +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                           breaks = seq(0, 1, 0.25)) +
        scale_fill_manual(values = cell_type_colors, name = "Cell Type",
                          drop = TRUE) +
        labs(x = NULL, y = "Proportion") +
        theme_bw(base_size = 18) +
        theme(strip.background   = element_blank(),
              axis.text.x        = element_blank(),
              axis.ticks.x       = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor   = element_blank(),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.margin = margin(5, 5, 0, 5)) +
        guides(fill = guide_legend(nrow = 2))

    ## Diagnosis annotation strip
    p_dx <- ggplot(donor_group_annot, aes(x = brnum, y = 1, fill = donor_group)) +
        geom_tile(width = 0.75) +
        scale_fill_manual(values = donor_group_colors, name = "Donor Group") +
        labs(x = "Sample") +
        theme_void(base_size = 18) +
        theme(axis.title.x = element_text(),
              plot.margin  = margin(0, 5, 5, 5))

    ## Combine with patchwork
    combined <- p_main / p_dx +
        plot_layout(heights = c(20, 1), guides = "collect") &
        theme(legend.position = "bottom",
              legend.box      = "horizontal")

    return(combined)
}

plot_box <- function(cell_prop){

    celltype_order <- c(
    "Astro", "Excit", "Immune", "Inhib", "Micro", 
    "D1-SPN", "D2-SPN", "Mural", "Oligo", "OPC"
    )

    cell_prop <- cell_prop |>
        mutate(cell_type = factor(cell_type, levels = celltype_order))

    donor_group_colors <- c(
    "BA" = "#a52a2a66",
    "WA" = "#0000ff66"
    )

    p <- ggplot(cell_prop, aes(x = cell_type, y = proportion)) +
        geom_boxplot(
          color= "black",
          fill = "grey60",
          width = 0.7,
          alpha = 0.8,
          outlier.shape = NA
        ) +
        geom_point(
          aes(color = donor_group),
          position = position_jitterdodge(
            jitter.width = 0.2,
            dodge.width = 0.4
          ),
          size = 1,
          alpha = 0.6
        ) +
        scale_color_manual(values = donor_group_colors, name = "Donor Group") +
        labs(x = "Cell Type", y = "Proportion") +
        theme_bw(base_size = 18) +
        theme(strip.background   = element_blank(),
              axis.text.x = element_text(angle = 35, hjust = 1),
              panel.grid.major.x = element_blank(),
              panel.grid.minor   = element_blank(),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.margin = margin(5, 5, 0, 5))
    
    return(p)
}

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}

## Main
out_path <- here("inputs/cell_proportions/_m")

if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

                                       # Generate phenotype data
pheno_file_path <- here("inputs/phenotypes/_m/phenotypes-all.tsv")

for (tissue in c("caudate", "dlpfc", "hippocampus")) {

                                        # Get sample list for valid ids
  samples_fn <- here("vmr-analysis/all_individuals/", paste0(tissue), "_m/samples.txt")
  samples <- fread(samples_fn, header = F, col.names = c("brnum", "FID")) %>%
    mutate(region = tissue)

                                        # Filter pheno file
  pheno_df <- clean_pheno(pheno_file_path, samples)

                                        # Get decovolution results
  cell_prop_fn <- here(paste0("inputs/cell_proportions/_m/music-proportions-", tissue, ".tsv"))
  cell_prop <- fread(cell_prop_fn, header = TRUE) %>%
    mutate(region = tissue)

                                        # Merge
  prop_df <- pheno_df %>%
    inner_join(cell_prop, by = c("brnum" = "sample_id"))

  summarize_prop(prop_df, tissue, out_path)

  p_stack <- plot_stack(prop_df)
  stacked_fn <- file.path(out_path, paste0("music_proportions_stacked-", tissue))
  save_plot(p_stack, stacked_fn, 12, 6)

  p_box <- plot_box(prop_df)
  box_fn <- file.path(out_path, paste0("music_proportions_box-", tissue))
  save_plot(p_box, box_fn, 8, 6)
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()