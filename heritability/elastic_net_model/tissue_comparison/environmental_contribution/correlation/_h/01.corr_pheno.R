#### Associate environmental phenotypes with low-heritability VMRs ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(tidyr)
})

## Function
filter_sites <- function(enet) {
  vmr <- na.omit(enet)
  vmr <- vmr %>% 
    select(chrom, start, end, h2_unscaled, r_squared_cv) %>% 
    mutate(h2_category = case_when(
      r_squared_cv <= 0.75 ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
    ),
    h2_category = factor(h2_category,
                         levels = c("Heritable", "Non-heritable", "Low prediction"))
    ) %>% 
    filter(h2_category == c("Non-heritable", "Low prediction"))
  
  return(vmr)
}

clean_pheno <- function(pheno_file_path, tissue, vars){
  pheno_df <- fread(pheno_file_path, header = TRUE) |>
    select(all_of(vars_to_include)) |>
    filter(agedeath > 17, region == tissue) |>
    mutate_if(is.character, as.factor)
  
  return(pheno_df)
}

merge_meth <- function(meth_files){
  meth_list <- vector("list", length(meth_files))
  for (i in seq_along(meth_files)) {
    file  <- meth_files[i]
                                        # Get pos from filename
    chr   <- basename(dirname(file))
    pos   <- strsplit(sub("_meth\\.phen$", "", basename(file)), "_")[[1]]
    df    <- fread(meth_files[i], select = c("V1", "V3"))
    colnames(df) <- c("brnum", "meth")
    
                                        # Add in vmr pos
    df <- df %>%
      mutate(chr   = as.integer(sub("chr_","", chr)),
             start = as.integer(pos[1]),
             end   = as.integer(pos[2]),
             feature_id = paste0("VMR", i))
    meth_list[[i]] <- df
  }
                                        # Bind meth matrix
  meth_df <- rbindlist(meth_list, use.names = TRUE, fill = TRUE)
  
  return(meth_df)
}

save_plot <- function(p, fn, w=8, h=6) {
  for (ext in c(".png", ".pdf")) {
    ggsave(filename = paste0(fn, ext), plot = p, width = w, height = h)
  }
}

## Main
tissue <- c("caudate")

out_path <- here("heritability", "elastic_net_model", "tissue_comparison",
                "environmental_contribution", "correlation", "_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

                                        # Read in summary table
enet_file <- here("heritability/elastic_net_model/", 
                  paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
enet <- read.table(enet_file, sep = "\t", header = TRUE)

vmr <- filter_sites(enet)

                                        # Define variables of interest
vars_to_include <- c(
  "brnum", "agedeath", "sex","primarydx","region","smoking","codeine","morphine",
  "cocaine","ethanol","antipsychotics","nicotine","amphetamines",
  "education","marital_status","hx_sexual_abuse","hx_physical_abuse",
  "hx_other_trauma","hx_military_service","fsiq"
)

                                        # Get phenotype data
pheno_file_path <- here("inputs/phenotypes/_m/phenotypes-AA.tsv")
pheno <- clean_pheno(pheno_file_path, tissue, vars_to_include)

                                        # Get meth matrix
meth_file_path <- here("heritability/caudate/_m/vmr")
meth_files     <- list.files(path = meth_file_path, pattern = "_meth\\.phen$", 
                         recursive = TRUE, full.names = TRUE)
meth_df        <- merge_meth(meth_files)

                                        # Merge with h2 groups and 
                                        # environmental variables
groups <- meth_df |>
  inner_join(vmr, "feature_id")
merged <- groups |>
  inner_join(pheno, "brnum") |>
  arrange(feature_id)

                                        # Write df to file
out_table <- file.path(out_path, paste0("vmr_env_assoc-AA.tsv.gz"))
fwrite(merged, out_table, sep = "\t", quote = FALSE)

                                        # Basic boxplot for exploratory
p <- ggplot(data=subset(merged, !is.na(smoking)), aes(x = smoking, y = meth)) +
  facet_wrap(~h2_category) +
  geom_boxplot()
p

merged %>% 
  group_by(smoking, h2_category) %>%
  summarize(
    count = n(),
    mean = mean(meth),
    sd = sd(meth)
  )

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()