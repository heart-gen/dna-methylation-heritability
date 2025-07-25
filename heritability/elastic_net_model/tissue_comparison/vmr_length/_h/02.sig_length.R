#### Check for length normality and significance ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(rstatix)
  library(ggpubr)
  library(ggplot2)
})

## Function 
filter_sites <- function(enet) {
  vmr <- na.omit(enet)
  vmr <- vmr %>%
    mutate(h2_category = case_when(
      r_squared_cv <= 0.75 ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
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

save_plot <- function(p, fn, w, h){
  for(ext in c('.png', '.pdf')){
    ggsave(file=paste0(fn, ext), plot=p, width=w, height=h)
  }
}

## Main
tissues <- c("caudate", "hippocampus", "dlpfc")

out_path <- here("heritability/elastic_net_model/tissue_comparison/vmr_length/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

vmr_all    <- list()

for (tissue in tissues) {
  # Read in summary table
  enet_file <- here("heritability/elastic_net_model/", 
                    paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
  enet <- read.table(enet_file, sep = "\t", header = TRUE)
  
  vmr <- filter_sites(enet)
  vmr <- cal_vmr_length(vmr)
  
  # Store vmrs across all tissues
  vmr$tissue <- tissue
  vmr_all[[tissue]] <- vmr
}

vmr_all <- bind_rows(vmr_all)
vmr_all$group <- interaction(vmr_all$tissue, vmr_all$h2_category, sep = "-")

heritability_test <- vmr_all %>%
  group_by(tissue) %>%
  pairwise_wilcox_test(log10_length ~ h2_category, p.adjust.method = "fdr") %>%
  add_y_position()

print(heritability_test)

heritability_colors <- c(
  "Heritable" = "#497C8A",
  "Non-heritable" = "#8CA77B",
  "Low prediction" = "#E3A27F"
)

p_heritability <- ggviolin(vmr_all, x = "h2_category",
                           y = "log10_length", fill = "h2_category", 
                           color = "h2_category", facet.by = "tissue", 
                           add = c("jitter", "median"), 
                           alpha = 0.5, palette = heritability_colors, 
                           trim = FALSE
) +
  stat_pvalue_manual(heritability_test, label = "p.adj.signif", tip.length = 0.01) +
  theme_pubr(base_size = 15, border = TRUE) +
  labs(
    title = "VMR Length Differences Across Heritability Categories",
    x = "Tissue", y = "log10(VMR Length)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

tissue_test <- vmr_all %>%
  group_by(h2_category) %>%
  pairwise_wilcox_test(log10_length ~ tissue, p.adjust.method = "fdr") %>%
  add_y_position()

print(tissue_test)

tissue_colors <- c(
  "caudate" = "#7372A6",
  "dlpfc" = "#B36F61",
  "hippocampus" = "#C5AC47"
)

p_tissue <- ggviolin(vmr_all, x = "tissue",
                     y = "log10_length", fill = "tissue", color = "tissue",
                     facet.by = "h2_category", add = c("jitter", "median"), 
                     alpha = 0.5, palette = tissue_colors, 
                     trim = FALSE
) +
  stat_pvalue_manual(tissue_test, label = "p.adj.signif", tip.length = 0.01) +
  theme_pubr(base_size = 15, border = TRUE) +
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
save_plot(p_tissue, fn_tissue, 14, 8)
save_plot(p_heritability, fn_heritability, 14, 8)

