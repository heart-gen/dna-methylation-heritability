## This script combines outputs from the elastic net model

## --- Main Script --- ##
                                        # Retrieve variables
ld_decay <- Sys.getenv("ld_decay")
method   <- Sys.getenv("METHOD", unset = "boosting_hybrid")

                                        # Function
read_data <- function(fn) {
    return(readr::read_table(fn, show_col_types=FALSE))
}

                                        # Loop through method-specific results directories
num_samples <- 200
for (dir_name in c("summary", "h2", "betas")) {
    outfile    <- paste("simulation", num_samples, ld_decay, method, dir_name,
                        "elastic-net.tsv", sep="_")
    dir_path   <- paste0(dir_name, "_", method)
    file_names <- list.files(dir_path, pattern="*.tsv$", full.names=TRUE)
    purrr::map_dfr(file_names, read_data) |>
        dplyr::mutate(PopSize = num_samples, LD_Decay = ld_decay,
                      Method = method,
                      ID = as.numeric(gsub("pheno_", "", pheno_id))) |>
        dplyr::arrange(ID) |> dplyr::select(-ID) |>
        data.table::fwrite(file=outfile, sep="\t")
}

## --- Reproducibility --- ##
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
