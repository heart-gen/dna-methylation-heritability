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

category_colors <- c(
  "Heritable" = "#497C8A",
  "Non-heritable" = "#8CA77B",
  "Low prediction" = "#E3A27F"
)

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

dunn <- vmr_all %>%
  dunn_test(log10_length ~ group, p.adjust.method = "fdr")

print(dunn)

# Generate y positions for labels
dunn_df <- dunn %>%
  mutate(y.position = max(vmr_all$log10_length, na.rm = TRUE) + seq(0.1, by = 0.1, length.out = n()))

# Create significance labels
dunn_df <- dunn_df %>%
  mutate(p.signif = cut(p.adj,
                        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                        labels = c("***", "**", "*", "ns")))

p <- ggplot(vmr_all, aes(x = tissue, y = log10_length, fill = h2_category)) +
  geom_boxplot(position = position_dodge(0.8), outlier.shape = NA) +
  geom_jitter(aes(color = h2_category), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
              size = 0.5, alpha = 0.3) +
  scale_fill_manual(values = category_colors) +
  scale_color_manual(values = category_colors) +
  stat_pvalue_manual(
    dunn_df,
    label = "p.signif",
    y.position = "y.position",
    xmin = "group1",
    xmax = "group2",
    tip.length = 0.01
  ) +
  theme_pubr(base_size = 14, border = TRUE) +
  labs(
    title = "VMR Length Differences Across Tissues and Heritability Categories",
    x = NULL,
    y = "log10(VMR Length)",
    fill = "Heritability Category",
    color = "Heritability Category"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

fn_box <- file.path(out_path, "VMR_length_comparisons")
save_plot(p, fn_box, 14, 8)