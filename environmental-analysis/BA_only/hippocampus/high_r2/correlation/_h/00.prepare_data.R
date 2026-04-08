#### Create matrix with methylation values and environmental variables ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

## Function
filter_sites <- function(enet) {
  vmr <- na.omit(enet)
  vmr <- vmr %>%
    dplyr::select(chrom, start, end, h2_unscaled, r_squared_cv) %>%
    mutate(h2_category = case_when(
      r_squared_cv <= 0.75 ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
    ),
    h2_category = factor(h2_category,
                         levels = c("Heritable", "Non-heritable", "Low prediction"))
    )

  return(vmr)
}

clean_pheno <- function(pheno_file_path, tissue, vars){
  pheno_df <- fread(pheno_file_path, header = TRUE) |>
    dplyr::select(all_of(vars)) |>
    filter(agedeath > 17, region == tissue) |>
    mutate(education = case_when(
      education %in% c("7th", "8th", "Less than 7th",
                       "9th", "10th", "11th", "12th") ~ "less_than_hs",
      education %in% c("H.S. diploma","GED") ~ "hs",
      education %in% c("1 yr college", "3 yrs college", "Associate's or 2 yrs college",
                       "Bachelor's", "Master's", "JD", "PhD") ~ "more_than_hs"
    )
    ) |>
    mutate(marital_status = case_when(
      marital_status %in% c("Single") ~ "single",
      marital_status %in% c("Married") ~ "married",
      marital_status %in% c("Divorced", "Separated", "Widowed") ~ "previously_married"
    )
    ) |>
    mutate(any_trauma_hx = dplyr::case_when(
      hx_sexual_abuse | hx_physical_abuse | hx_other_trauma | hx_military_service ~ 1L,
      is.na(hx_sexual_abuse) & is.na(hx_physical_abuse) &
        is.na(hx_other_trauma) & is.na(hx_military_service) ~ NA_integer_,
      TRUE ~ 0L
    )) |>
    rename(age = agedeath, dx = primarydx) |>
    mutate_if(is.character, as.factor)

    print(colnames(pheno_df))

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

## Main
tissue <- c("hippocampus")

out_path <- here("environmental-analysis", "BA_only",
                 paste0(tissue, "/high_r2/correlation/_m"))
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

                                        # Read in summary table
enet_file <- here("heritability/elastic_net_model/BA_only/",
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

                                        # Add global ances
f_ances   <- here("inputs", "genetic-ancestry",
                  "structure.out_ancestry_proportion_raceDemo_compare")
ances     <- fread(f_ances) %>% filter(group == "AA")
pheno <- left_join(pheno, ances, by = c("brnum" = "id")) %>%
  rename(afr_ances = Afr)

                                        # Get meth matrix
meth_file_path <- here("vmr-analysis/hippocampus/_m/vmr")
meth_files     <- list.files(path = meth_file_path, pattern = "_meth\\.phen$",
                         recursive = TRUE, full.names = TRUE)
meth_df        <- merge_meth(meth_files)
meth_df <- meth_df %>% mutate(chr = as.character(chr))

                                        # Merge with h2 groups and
                                        # environmental variables
groups <- meth_df |>
  inner_join(vmr, by = c("chr" = "chrom", "start", "end"))
merged <- groups |>
  inner_join(pheno, "brnum") |>
  arrange(feature_id)

                                        # Write df to file
out_table <- file.path(out_path, paste0("vmr_env_assoc-AA.tsv.gz"))
fwrite(merged, out_table, sep = "\t", quote = FALSE)

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()