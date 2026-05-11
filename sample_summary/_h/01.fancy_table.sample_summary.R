## Summarize sample data

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(gtsummary)
    library(data.table)
})

options(gt.html_tag_check = FALSE)

## Function
save_table <- function(pp, fn){
    for(ext in c(".tex", ".rtf")){
        gt::gtsave(as_gt(pp), filename=paste0(fn,ext))
    }
}

clean_pheno <- function(pheno_file_path, samples_to_include){
  pheno_df <- fread(pheno_file_path, header = TRUE) |>
    inner_join(samples_to_include, by = c("brnum", "region")) |>
    filter(agedeath > 17) |>
    mutate(region = recode(region,
                           "caudate" = "Caudate",
                           "dlpfc" = "DLPFC",
                           "hippocampus" = "Hippocampus")) |>
    mutate_if(is.character, as.factor)
    
  return(pheno_df)
}

## Main
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

samples_to_include <- bind_rows(valid_samples)

pheno_df <- clean_pheno(pheno_file_path, samples_to_include)

data.table::fwrite(pheno_df, file="phenotype_data.tsv", sep="\t",
                   row.names=FALSE, col.names=TRUE)

                                         # Generate pretty tables
fn_AA <- "sample_breakdown_AA.table"
pp_AA <- pheno_df |>
    filter(race == "AA") |>
    select(agedeath, sex, primarydx, region) |>
    mutate(sex = factor(sex, labels=c("Female", "Male")),
           primarydx = factor(forcats::fct_drop(primarydx), labels=c("CTL", "SCZ"))) |>
    tbl_summary(by="region", missing="no",
                label = list(agedeath ~ "Age", sex ~ "Sex", primarydx ~ "Dx"),
                statistic=all_continuous() ~ c("{mean} ({sd})")) |>
    modify_header(all_stat_cols() ~ "**{level}**<br>N = {n}") |>
    modify_spanning_header(all_stat_cols() ~ "**Brain Region**") |>
    bold_labels() |> italicize_levels()
print(pp_AA)
save_table(pp_AA, fn_AA)

fn_all <- "sample_breakdown_all.table"
pp_all <- pheno_df |>
  select(agedeath, race, sex, primarydx, region) |>
  mutate(sex = factor(sex, labels=c("Female", "Male")),
         primarydx = factor(forcats::fct_drop(primarydx), labels=c("CTL", "SCZ"))) |>
  tbl_strata(
    strata=region,
    ~.x %>% 
      tbl_summary(by="race", missing="no",
                  label = list(agedeath ~ "Age", sex ~ "Sex", primarydx ~ "Dx"),
                  statistic=all_continuous() ~ c("{mean} ({sd})")) |>
      modify_header(all_stat_cols() ~ "**{level}**<br>N = {n}") |>
      modify_spanning_header(all_stat_cols() ~ "**Brain Region**") |>
      bold_labels() |> italicize_levels()
  )
print(pp_all)
save_table(pp_all, fn_all)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
