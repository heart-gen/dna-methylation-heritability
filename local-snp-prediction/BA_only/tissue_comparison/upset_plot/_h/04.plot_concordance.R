## This script plots concordance for VMR sharing across brain regions

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

## Function
count_intersections <- function(fn){
  vmrs <- fread(fn)
  return(nrow(vmrs))
}

## Main
out_path <- here("heritability/elastic_net_model/BA_only/tissue_comparison/upset_plot/_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

option_flags <- c("F_0.25", "F_0.5", "F_0.75", "f_0.25", "f_0.5", "f_0.75")
h2_categories <- c("heritable", "non-heritable", "low_prediction", "all")

results <- list()

for (flag in option_flags) {
  for (h2_cat in h2_categories){
    
    # Get overlapping VMRs
    sets <- c(
      Caudate = count_intersections(paste0("./", flag, "/sets/caudate_specific_", h2_cat, ".bed")),
      DLPFC = count_intersections(paste0("./", flag, "/sets/dlpfc_specific_", h2_cat, ".bed")),
      Hippocampus = count_intersections(paste0("./", flag, "/sets/hippocampus_specific_", h2_cat, ".bed")),
      `Caudate&Hippocampus` = count_intersections(paste0("./", flag, "/sets/caudate_hippocampus_overlap_", h2_cat, ".bed")),
      `Caudate&DLPFC` = count_intersections(paste0("./", flag, "/sets/caudate_dlpfc_overlap_", h2_cat, ".bed")),
      `Hippocampus&DLPFC` = count_intersections(paste0("./", flag, "/sets/hippocampus_dlpfc_overlap_", h2_cat, ".bed")),
      `Caudate&Hippocampus&DLPFC` = count_intersections(paste0("./", flag, "/sets/3tissues_overlap_", h2_cat, ".bed.tmp"))
    )
    
    df <- data.frame(
      flag = flag,
      h2_category = h2_cat,
      set = names(sets),
      count = as.numeric(sets)
    )
    
    results[[paste(flag, h2_cat, sep = "_")]] <- df
    
  }
}

df <- bind_rows(results)
fwrite(df, file.path(out_path, "overlap_summary.tsv"), sep = "\t")

# Get proportions 
df_prop <- df %>%
  group_by(flag, h2_category) %>% 
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Plot proportions 
p <- ggplot(df_prop, aes(x = set, y = prop, fill = h2_category)) +
  geom_bar(
    stat = "identity",
    position = position_stack(), 
    color = "black"
  ) +
  scale_fill_manual(
    values = c("all" = "grey",
               "heritable" = "#497C8A",
               "non-heritable" = "#8CA77B",
               "low_prediction" = "#E3A27F")
  ) +
  labs(
    x = "Overlap Set",
    y = "Proportion",
    fill = "Category"
  ) +
  facet_wrap(~flag) +
  theme_minimal(base_size = 20) +
  #font("xy.title", face = "bold", size = 14) + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.title = element_text(hjust = 0.5),
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()