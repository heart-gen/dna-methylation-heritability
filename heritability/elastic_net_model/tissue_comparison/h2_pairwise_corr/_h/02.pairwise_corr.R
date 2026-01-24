#### Correlate heritability estimates across brain regions ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(GenomicRanges)
  library(data.table)
  library(ggpubr)
})

## Function
save_plot <- function(p, fn, w, h){
  for(ext in c(".pdf", ".png")){
    ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
  }
}

get_overlaps <- function(vmr_all, h2_cat, tissue1, tissue2){
  overlap_fn <- paste0("F_0.5/", tissue1, "_", tissue2, "_overlap_0.5.bed")
  shared_tissue1 <- fread(overlap_fn, select = 1:3, 
                          col.names = c("chr", "start", "end")) |>
                    mutate(VMR_id = row_number())
  shared_tissue2 <- fread(overlap_fn, select = 4:6, 
                          col.names = c("chr", "start", "end")) |>
                    mutate(VMR_id = row_number())
  
  vmr_tissue1 <- vmr_all %>%
    filter(tissue == tissue1, h2_category == h2_cat)
  
  vmr_tissue2 <- vmr_all %>%
    filter(tissue == tissue2, h2_category == h2_cat)
  
  merged_tissue1 <- vmr_tissue1 %>%
    right_join(shared_tissue1, by = c("chrom" = "chr", "start", "end")) %>%
    arrange(VMR_id) %>% na.omit()
  
  merged_tissue2 <- vmr_tissue2 %>%
    right_join(shared_tissue2, by = c("chrom" = "chr", "start", "end")) %>%
    arrange(VMR_id) %>% na.omit()
  
  merged <- merged_tissue1 %>%
    inner_join(merged_tissue2, by = "VMR_id")
  
  h2_df <- data.frame(
    h2_tissue1 = merged$h2_unscaled.x,
    h2_tissue2 = merged$h2_unscaled.y
  )
  
  return(h2_df)
}

spearman_corr <- function(h2_df, h2_cat, tissue1, tissue2, out_path){
  spearman <- h2_df %>% 
    summarise(
      spearman_rho = cor.test(h2_tissue1, h2_tissue2, method = "spearman")$estimate,
      spearman_p_value = cor.test(h2_tissue1, h2_tissue2, method = "spearman")$p.value,
      n = n()
    )
  print(spearman)
  write.csv(spearman, 
            file = file.path(out_path, 
                             paste0(gsub(" ", "_", tolower(h2_cat)), "_h2_corr_", 
                                    tissue1, "_", tissue2, ".csv")), 
            row.names = FALSE)
}

plot_corr <- function(h2_df, h2_cat, tissue1, tissue2, output_path){
  xlab = paste0("h2\n", tissue1)
  ylab = paste0("h2\n", tissue2)
  
  fn = file.path(output_path, paste(gsub(" ", "_", tolower(h2_cat)), "h2_corr",
                                    tissue1, tissue2, sep="_"))
  
  pp <- ggscatter(h2_df, x = "h2_tissue1", y = "h2_tissue2", add = "reg.line", 
                  size = 1, xlab = xlab, ylab = ylab, panel.labs.font=list(face="bold"),
                  add.params=list(color = "blue", fill = "lightgray"),
                  conf.int = TRUE, cor.coef = TRUE, cor.coef.size = 3,
                  cor.method = "spearman", cor.coeff.args=list(label.sep="\n"), ncol = 4) +
    theme_pubr(base_size=18) + 
    labs(title = paste("h2 correlation:", h2_cat))
  save_plot(pp, fn, 6, 6)
}

## Main
tissues <- c("caudate", "hippocampus", "dlpfc")

out_path <- here("heritability/elastic_net_model/tissue_comparison/h2_pairwise_corr/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

vmr_all    <- list()

for (tissue in tissues) {
  # Read in summary table
  enet_file <- here("heritability/elastic_net_model/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
  vmr <- read.table(enet_file, sep = "\t", header = TRUE, 
                    colClasses = c(chrom = "character")) |> na.omit()
  
  # Store vmrs across all tissues
  vmr$tissue <- tissue
  vmr_all[[tissue]] <- vmr
}

vmr_all <- bind_rows(vmr_all)

vmr_all <- vmr_all %>%
  mutate(h2_category = case_when(
    r_squared_cv <= 0.75 ~ "Low prediction",
    h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
    h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
  ),
  h2_category = factor(h2_category, levels = c("Heritable", 
                                               "Non-heritable", 
                                               "Low prediction"))
  )

for (h2_cat in c("Heritable", "Non-heritable", "Low prediction")){
  for (tissue1 in c("caudate", "hippocampus")){
    for (tissue2 in c("hippocampus", "dlpfc")){
      if (tissue1 != tissue2){
        h2_df <- get_overlaps(vmr_all, h2_cat, tissue1, tissue2)
        spearman_corr(h2_df, h2_cat, tissue1, tissue2, out_path)
        plot_corr(h2_df, h2_cat, tissue1, tissue2, out_path)
      }
    }
  }
}

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()