#### Check for length normality and significance ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(rstatix)
  library(ggpubr)
  library(ggplot2)
  library(lme4)
  library(emmeans)
  library(data.table)
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

cal_vmr_length <- function(vmr) {
  vmr <- vmr %>%
    mutate(length = end - start,
           log10_length = log10(length))
  return(vmr)
}

get_samples <- function(tissue) {
  sample_fn <- here("vmr-analysis/", paste0(tissue, "/_m/samples.txt"))
  samples <- fread(sample_fn, header = FALSE)
  n_samples <- nrow(samples)
  return(n_samples)
}

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}

## Main
tissues <- c("caudate", "hippocampus", "dlpfc")

out_path <- here("heritability/elastic_net_model/BA_only/tissue_comparison/vmr_length/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

vmr_all    <- list()

for (tissue in tissues) {
  # Read in summary table
  enet_file <- here("heritability/elastic_net_model/BA_only/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
  enet <- read.table(enet_file, sep = "\t", header = TRUE)
  
  vmr <- filter_sites(enet)
  vmr <- cal_vmr_length(vmr)
  
  # Store vmrs across all tissues
  vmr$tissue <- tissue
  vmr_all[[tissue]] <- vmr
}

vmr_all <- bind_rows(vmr_all)

# Test significance with mixed linear model
contrasts_h2 <- list()

for(t in unique(vmr_all$tissue)) {
  
  n_samples <- get_samples(t)
  df <- vmr_all %>% filter(tissue == t) %>%
    mutate(n_samples = n_samples)
  
  # Fit mixed model per tissue
  lmm <- lmer(log10_length ~ h2_category + n_samples + (1 | chrom), data = df)
  y_pos = max(df$log10_length, na.rm = TRUE) + 0.6
  
  # Get pairwise contrasts
  contrast_df <- emmeans(lmm, ~ h2_category) %>%
    contrast(method = "pairwise", adjust = "bonferroni") %>%
    as.data.frame() %>%
    mutate(
      tissue = t,
      group1 = gsub("[()]", "", sapply(strsplit(contrast, " - "), `[`, 1)),
      group2 = gsub("[()]", "", sapply(strsplit(contrast, " - "), `[`, 2)),
      y.position = y_pos,
      p.adj.signif = case_when(
        p.value <= 0.001 ~ "***",
        p.value <= 0.01  ~ "**",
        p.value <= 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
  
  contrasts_h2[[t]] <- contrast_df
}

heritability_test <- bind_rows(contrasts_h2)
write.csv(heritability_test, 
            file = file.path(out_path, 
                             paste0("vmr_length_heritability_comparisons.csv")), 
            row.names = FALSE)

heritability_colors <- c(
  "Heritable" = "#497C8A",
  "Non-heritable" = "#8CA77B",
  "Low prediction" = "#E3A27F"
)

p_heritability <- ggviolin(vmr_all, x = "h2_category",
                           y = "log10_length", fill = "h2_category", 
                           facet.by = "tissue", add = "boxplot", 
                           alpha = 0.5, palette = heritability_colors, 
                           trim = FALSE
) +
  stat_pvalue_manual(heritability_test, label = "p.adj.signif", tip.length = 0.01) +
  theme_pubr(base_size = 18, border = TRUE) +
  labs(
    title = "VMR Length Differences Across Heritability Categories",
    x = "Tissue", y = "log10(VMR Length)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

contrasts_tissue <- list()

for(h in unique(vmr_all$h2_category)) {
  
  n_samples <- get_samples(t)
  df <- vmr_all %>% filter(h2_category == h) %>%
    mutate(n_samples = n_samples)
  
  # Fit mixed model per h2 category 
  lmm <- lmer(log10_length ~ tissue + n_samples + (1 | chrom), data = df)
  y_pos = max(df$log10_length, na.rm = TRUE) + 0.6
  
  # Get pairwise contrasts
  contrast_df <- emmeans(lmm, ~ tissue) %>%
    contrast(method = "pairwise", adjust = "bonferroni") %>%
    as.data.frame() %>%
    mutate(
      h2_category = h,
      group1 = gsub("[()]", "", sapply(strsplit(contrast, " - "), `[`, 1)),
      group2 = gsub("[()]", "", sapply(strsplit(contrast, " - "), `[`, 2)),
      y.position = y_pos,
      p.adj.signif = case_when(
        p.value <= 0.001 ~ "***",
        p.value <= 0.01  ~ "**",
        p.value <= 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
  
  contrasts_tissue[[h]] <- contrast_df
}

tissue_test <- bind_rows(contrasts_tissue)
write.csv(tissue_test, 
            file = file.path(out_path, 
                             paste0("vmr_length_tissue_comparisons.csv")), 
            row.names = FALSE)

tissue_colors <- c(
  "caudate" = "#7372A6",
  "dlpfc" = "#B36F61",
  "hippocampus" = "#C5AC47"
)

p_tissue <- ggviolin(vmr_all, x = "tissue",
                     y = "log10_length", fill = "tissue",
                     facet.by = "h2_category", add = "boxplot", 
                     alpha = 0.5, palette = tissue_colors, 
                     trim = FALSE
) +
  stat_pvalue_manual(tissue_test, label = "p.adj.signif", tip.length = 0.01) +
  theme_pubr(base_size = 18, border = TRUE) +
  labs(
    title = "VMR Length Differences Across Tissues",
    x = "Tissue", y = "log10(VMR Length)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

fn_heritability <- file.path(out_path, "VMR_length_heritability_comparisons")
fn_tissue       <- file.path(out_path, "VMR_length_tissue_comparisons")
save_plot(p_tissue, fn_tissue, 10, 7)
save_plot(p_heritability, fn_heritability, 10, 7)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

