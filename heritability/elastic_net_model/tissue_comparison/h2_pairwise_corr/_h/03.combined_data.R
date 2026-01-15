## This script combines outputs from correlation analysis

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(stringr)
})

## Function
read_data <- function(fn) {
    vars <- str_match(basename(fn), "(.*)_h2_corr_(.*)_(.*)\\.csv")
    
    df <- fread(fn) |> 
      mutate(h2_category = vars[2],
             tissue1 = vars[3],
             tissue2 = vars[4])
    
    return(df)
}

## Main
                                        # Loop through results directory

outfile    <- "vmr_h2_corr_spearman_summary.tsv"
file_names <- list.files(pattern = "*.csv$", full.names = TRUE)
purrr::map_dfr(file_names, read_data) |>
    fwrite(file=outfile, sep="\t")

## --- Reproducibility --- ##
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
