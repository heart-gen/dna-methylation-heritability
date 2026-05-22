#### Plot SNP PC1 and PC2 for all individuals ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(ggpubr)
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
            race = recode(race, "AA" = "BA", "EA" = "WA")) |>
    mutate_if(is.character, as.factor)
    
  return(pheno_df)
}

plot_pc <- function(pheno_df){

  # Define palette
  donor_group_colors <- c(
    "BA" = "#a52a2a66",
    "WA" = "#0000ff66"
  )
  
  p <- ggscatter(pheno_df, x = "snpPC1", y = "snpPC2",
                      size = 1.5, alpha = 0.75,
                      xlab = "SNP PC1", ylab = "SNP PC2", 
                      color = "race",
                      ggtheme = theme_pubr(base_size = 16, border = TRUE)
  ) +
    facet_wrap(~region, scales = "free_y") +
    scale_color_manual(values = donor_group_colors) +
    labs(color = NULL) +
    font("xy.title", face = "bold", size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )

  return(p)
}

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}

## Main
out_path <- here("qc_analysis/replication_cohort/_m")

if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}
                                        # Generate phenotype data
pheno_file_path <- here("inputs/phenotypes/_m/phenotypes-all.tsv")

                                        # Get sample list for valid ids
valid_samples <- list()

for (tissue in c("caudate", "dlpfc", "hippocampus")) {
  samples_fn <- here("vmr-analysis/all_individuals/", paste0(tissue), "_m/samples.txt")
  samples <- fread(samples_fn, header = F, col.names = c("brnum", "FID")) %>%
    mutate(region = tissue)
  
  valid_samples[[tissue]] <- samples
}

                                        # Combine across brain regions
samples_to_include <- bind_rows(valid_samples)

                                        # Filter pheno file
pheno_df <- clean_pheno(pheno_file_path, samples_to_include)

                                        # Plot PC1 vs PC2
p <- plot_pc(pheno_df)

                                        # Write to file
fn_pc <- file.path(out_path, "SNP_PC_AA_EA")
save_plot(p, fn_pc, 9, 4.5)


#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

