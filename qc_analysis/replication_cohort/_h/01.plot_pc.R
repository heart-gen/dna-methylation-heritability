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
            pop = recode(race, "AA" = "BA", "EA" = "WA")) |>
    mutate_if(is.character, as.factor)
    
  return(pheno_df)
}

plot_pc <- function(pheno_df, ref_df){

  # Define palette
  donor_group_colors <- c(
    "BA" = "#a52a2a66",
    "WA" = "#0000ff66"
  )

  # Define palette
  ref_group_colors <- c(
    "AFR" = "#e41a1c",
    "EUR" = "#377eb8",
    "EAS" = "#4daf4a",
    "SAS" = "#984ea3",
    "AMR" = "#ff7f00"
  )
  
  ggplot() +
    geom_point(
      data = ref_df,
      aes(x = snpPC1, y = snpPC2, shape = 17, color = super_pop),
      size = 1,
      alpha = 0.5
    ) +
    geom_point(
      data = pheno_df,
      aes(x = snpPC1, y = snpPC2, color = pop),
      size = 1.5,
      alpha = 0.8
    ) +
    facet_wrap(~region, scales = "free_y") +
    ggtheme = theme_pubr(base_size = 16, border = TRUE) +
    scale_color_manual(values = c(ref_group_colors, donor_group_colors)) +
    labs(color = NULL, x = "SNP PC1", y = "SNP PC2") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
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

                                        # Get 1000GP ref samples
ref_samples_fn <- "/projects/b1213/resources/1kGP/data_raw/integrated_call_samples_v3.20130502.ALL.panel"
ref_samples <- fread(ref_samples_fn, header = TRUE) |>
    col.names(c("IID", "pop", "super_pop", "gender"))

ref_pc_fn <- "/projects/b1213/resources/1kGP/GRCh38_phased_vcf/_m/1kGP.eigenvec"
ref_pc <- fread(ref_pc_fn, header = FALSE) |>
    col.names(c("FID", "IID", paste0("snpPC", 1:10))) |>
    mutate(dataset = "1000 Genomes") |>
    left_join(ref_samples, by = "IID")

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
p <- plot_pc(pheno_df, ref_pc)

                                        # Write to file
fn_pc <- file.path(out_path, "SNP_PC_AA_EA")
save_plot(p, fn_pc, 10, 8)


#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

