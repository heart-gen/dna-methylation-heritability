## Prepare covariates for FastQTL / tensorQTL analysis
## Requires the 'sva' package
suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(SummarizedExperiment)
    library('data.table')
})

#### Functions
format_pheno <- function(pheno_file_path) {
    pheno_df <- fread(pheno_file_path, header = TRUE) |>
        filter(race == "AA", agedeath >= 17, region == "caudate") |>
        select(brnum, agedeath, sex, primarydx)
    
                                        # One hot encode categorical covs
    one_hot_df <- model.matrix(~ primarydx + sex - 1, data = pheno_df) |>
        as.data.frame() |>
        select(-primarydxControl)
    
    pheno <- cbind(BrNum = pheno_df$brnum, 
                   one_hot_df, 
                   Age = pheno_df$agedeath)
  return(pheno)
}

get_covs <- function(pheno, sample_id) {
                                      # Subset samples
    sample_df <- data.table::fread(sample_id)
  
                                      # Extract covariates
    covs_df <- pheno |>
      inner_join(sample_df, by = "BrNum") |>
      tibble::column_to_rownames("BrNum")
    
                                      # Reformat covariates
    covs <- t(covs_df) |>
      as.data.frame() |>
      tibble::rownames_to_column("ID") |>
      select(c("ID", sample_df$BrNum))
    
    return(covs)
}

#### Main

                                       # Create output dir
output_path <- here("meqtl-analysis", "vmrs", "caudate", "_m")

if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

                                        # Get pheno file for covariates
pheno_file_path <- here("inputs", "phenotypes", "_m", "phenotypes-AA.tsv")
pheno <- format_pheno(pheno_file_path)


                                        # Extract covariates
sample_id <- here("meqtl-analysis", "vmrs", "caudate", "_m", 
                  "sample_id_to_brnum.tsv")
covs <- get_covs(pheno, sample_id)
out_covs   <- file.path(output_path, "vmrs.combined_covariates.txt")

                                        # Save file
data.table::fwrite(covs, out_covs, sep='\t')

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
