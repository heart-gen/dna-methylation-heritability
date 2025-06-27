## This script extracts phenotype information from the
## BSseq objects (DNAm).

load_methyl <- function(fn) {
    load(here::here("inputs/wgbs-data", fn))
    tissue      <- dirname(fn)
    sample_data <- bsseq::pData(BSobj) |>
        as.data.frame() |> dplyr::filter(race %in% c("AA")) |>
        dplyr::mutate_if(is.character, as.factor) |>
        dplyr::mutate(region = tissue)
    return(sample_data)
}

#### MAIN
                                        # Load methylation
fnames  <- c("caudate/Caudate_chr21_BSobj.rda",
             "hippocampus/Hippocampus_chr21BSobj_Genotypes.rda",
             "dlpfc/dlpfc_chr21BSobj_GenotypesDosage.rda")

purrr::map_dfr(fnames, load_methyl) |>
    dplyr::select(c("brnum", "agedeath", "sex", "race", "primarydx",
                    "pmi", "batch", "region")) |>
    data.table::fwrite("phenotypes-DNAm.tsv", sep="\t")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
