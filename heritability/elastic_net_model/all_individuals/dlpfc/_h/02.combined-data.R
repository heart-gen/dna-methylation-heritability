## This script combines outputs from the elastic net model

## --- Main Script --- ##
                                        # Retrieve variables
region  <- Sys.getenv("region")

                                        # Function
read_data <- function(fn) {
    return(readr::read_table(fn, show_col_types=FALSE) |> 
             dplyr::mutate(chrom = as.character(chrom))
           )
}

                                        # Loop through results directory
for (dir_name in c("summary", "h2", "betas")) {
    for (pop in c("AA", "EA")) {
        outfile    <- paste0(tolower(region), "_", dir_name, 
                            "_elastic-net_", pop, ".tsv")
        file_names <- list.files(dir_name,
                                 pattern=paste0(pop, ".*\\.tsv$"),
                                 full.names=TRUE)
        purrr::map_dfr(file_names, read_data) |>
    	dplyr::mutate(region = region) |>
        data.table::fwrite(file=outfile, sep="\t")
    }
}

## --- Reproducibility --- ##
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
