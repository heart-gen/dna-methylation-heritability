## This script combines outputs from the elastic net model

## --- Main Script --- ##
                                        # Retrieve variables
region  <- Sys.getenv("region")

                                        # Function
read_data <- function(fn) {
    return(readr::read_table(fn, show_col_types=FALSE))
}

                                        # Loop through results directory
for (dir_name in c("summary", "h2", "betas")) {
    outfile    <- paste(tolower(region), dir_name, "elastic-net.tsv", sep="_")
    file_names <- list.files(dir_name,pattern="*.tsv$",full.names=TRUE)
    purrr::map_dfr(file_names, read_data) |>
        dplyr::mutate(region = region) |>
        data.table::fwrite(file=outfile, sep="\t")
}

## --- Reproducibility --- ##
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
