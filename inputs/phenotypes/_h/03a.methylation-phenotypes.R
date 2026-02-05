## This script extracts phenotype information from the
## BSseq objects (DNAm).

load_methyl_rda <- function(fn) {
    load(here::here("inputs/wgbs-data", fn))
    tissue      <- dirname(fn)
    sample_data <- bsseq::pData(BSobj) |>
        as.data.frame() |> dplyr::filter(race %in% c("AA", "CAUC")) |>
        dplyr::mutate_if(is.character, as.factor) |>
        dplyr::mutate(region = tissue)
    sample_data$race <- gsub("CAUC", "EA", sample_data$race)
    return(sample_data)
}

load_methyl_h5 <- function(fn) {
    BSobj       <- HDF5Array::loadHDF5SummarizedExperiment(dir = fn)
    tissue      <- dirname(fn)
    sample_data <- bsseq::pData(BSobj) |>
        as.data.frame() |> dplyr::filter(race %in% c("AA", "CAUC")) |>
        dplyr::mutate_if(is.character, as.factor) |>
        dplyr::mutate(region = tissue)
    sample_data$race <- gsub("CAUC", "EA", sample_data$race)
    return(sample_data)
}

#### MAIN
                                        # Load methylation
fn1    <- "caudate/_m/caudate_chr21_BSobj.rda"
fnames <- c("hippocampus/_m/hippocampus_chr21_BSobj",
            "dlpfc/_m/dlpfc_chr21_BSobj")
df1 <- load_methyl_rda(fn1) |>
    dplyr::select(c("brnum", "agedeath", "sex", "race", "primarydx",
                    "pmi", "region"))
df2 <- purrr::map_dfr(fnames, load_methyl_h5) |>
    dplyr::select(c("brnum", "agedeath", "sex", "race", "primarydx",
                    "pmi", "region"))

dplyr::bind_rows(df1, df2) |>
    data.table::fwrite("phenotypes-DNAm-all.tsv", sep="\t")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
