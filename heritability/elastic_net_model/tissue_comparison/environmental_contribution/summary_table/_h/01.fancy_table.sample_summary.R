## Summarize environmental factors by brain region

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

clean_table <- function(pheno_file_path, vars){
    pheno_df <- fread(pheno_file_path, header = TRUE) |>
    select(all_of(vars_to_include)) |>
    filter(agedeath > 17) |>
    mutate(region = recode(region,
                           "caudate" = "Caudate",
                           "dlpfc" = "DLPFC",
                           "hippocampus" = "Hippocampus")) |>
  mutate(education_group = case_when(
      education %in% c("7th", "8th", "Less than 7th") ~ "Less than high school",
      education %in% c("9th", "10th", "11th", "12th") ~ "Some high school",
      education %in% c("H.S. diploma","GED") ~ "High School",
      education %in% c("1 yr college", "3 yrs college", "Associate's or 2 yrs college") ~ "Some college",
      education == "Bachelor's" ~ "Bachelor's",
      education == "Master's" ~ "Master's",
      education %in% c("JD", "PhD") ~ "Doctorate"
    )
    ) |>
    mutate_if(is.character, as.factor)

    return(pheno_df)
}

## Main
                                        # Get variables of interest
vars_to_include <- c(
  "agedeath", "sex","primarydx","region","smoking","codeine","morphine",
  "cocaine","ethanol","antipsychotics","nicotine","amphetamines",
  "education","marital_status","hx_sexual_abuse","hx_physical_abuse",
  "hx_other_trauma","hx_military_service","fsiq"
)

                                        # Generate phenotype data
pheno_file_path <- here("inputs/phenotypes/_m/phenotypes-AA.tsv")
pheno_df <- clean_table(pheno_file_path, vars_to_include)

                                         # Generate pretty tables
fn <- "environmental_factors.table"
pp <- pheno_df |> select(!c("agedeath", "education")) |>
    mutate(sex=factor(sex, labels=c("Female", "Male")),
           primarydx=factor(forcats::fct_drop(primarydx), labels=c("CTL", "SCZ"))) |>
    tbl_summary(by="region", missing="no",
                label = list(
                  sex ~ "Sex", primarydx ~ "Diagnosis", smoking ~ "Smoking",
                  codeine ~ "Codeine", morphine ~ "Morphine", cocaine ~ "Cocaine",
                  ethanol ~ "Ethanol", amphetamines ~ "Amphetamines",
                  antipsychotics ~ "Antipsychotics", nicotine ~ "Nicotine",
                  hx_military_service ~ "Military Service", hx_other_trauma ~ "Other Trauma",
                  hx_physical_abuse ~ "Physical Abuse", hx_sexual_abuse ~ "Sexual Abuse",
                  education_group ~ "Education", marital_status ~ "Marital Status",
                  fsiq ~ "FSIQ"),
                statistic=all_continuous() ~ c("{mean} ({sd})")) |>
    modify_header(all_stat_cols()~"**{level}**<br>N = {n}") |>
    modify_spanning_header(all_stat_cols()~"**Brain Region**") |>
    bold_labels() |> italicize_levels()
save_table(pp, fn)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
