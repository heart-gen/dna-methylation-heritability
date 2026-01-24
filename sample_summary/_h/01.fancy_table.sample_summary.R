## Summarize sample data

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(gtsummary)
})

options(gt.html_tag_check = FALSE)

## Function
save_table <- function(pp, fn){
    for(ext in c(".tex", ".rtf")){
        gt::gtsave(as_gt(pp), filename=paste0(fn,ext))
    }
}

clean_pheno <- function(pheno_file_path){
  pheno_df <- fread(pheno_file_path, header = TRUE) |>
    select(agedeath, sex, primarydx, region) |>
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
pheno_file_path <- here("inputs/phenotypes/_m/phenotypes-AA.tsv")
pheno_df <- clean_pheno(pheno_file_path)

data.table::fwrite(pheno_df, file="phenotype_data.tsv", sep="\t",
                   row.names=FALSE, col.names=TRUE)

                                         # Generate pretty tables
fn <- "sample_breakdown.table"
pp <- pheno_df |>
    mutate(sex = factor(sex, labels=c("Female", "Male")),
           primarydx = factor(forcats::fct_drop(primarydx), labels=c("CTL", "SCZ"))) |>
    tbl_summary(by="region", missing="no",
                label = list(agedeath ~ "Age", sex ~ "Sex", primarydx ~ "Dx"),
                statistic=all_continuous() ~ c("{mean} ({sd})")) |>
    modify_header(all_stat_cols() ~ "**{level}**<br>N = {n}") |>
    modify_spanning_header(all_stat_cols() ~ "**Brain Region**") |>
    bold_labels() |> italicize_levels()
print(pp)
save_table(pp, fn)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
