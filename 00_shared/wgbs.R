#### WGBS BSseq object loading (v2 revision) ####
##
## Loading a BSobj is not a plain load(), because the region trees do not agree
## on what a BSobj is (V14):
##
##   caudate      inputs/wgbs-data/caudate/_m/*.rda -- a DelayedArray view whose
##                baked-in HDF5 path is valid. M, Cov and coef all resolve.
##   dlpfc        the .rda names a backing file (dlpfc_assays.h5) that no longer
##                exists -- the owner renamed it to assays.h5 -- and its `coef`
##                assay points into a /tmp scratch dir from the smoothing
##                session, so the smoothed fits are unrecoverable from the .rda.
##   hippocampus  the .rda's M/Cov resolve, but `coef` is /tmp-backed and gone.
##
## For dlpfc and hippocampus there is a self-contained
## saveHDF5SummarizedExperiment() directory alongside the .rda that carries all
## three assays intact. v2 loads that in preference to the .rda. Verified
## 2026-08-15 for chr22 in both regions: rowRanges, colnames and colData$brnum
## are identical() to the .rda, and the directory form additionally has usable
## `coef` (hasBeenSmoothed TRUE).
##
## This matters because 00_prepare.R uses getMeth(type = "smooth"), which reads
## `coef`. An object missing it fails only once the pipeline reaches that call.
##
## AGENTS.md 9: Quest paths belong in configuration, so the per-region source
## lives in config/paths.yml, not in this file.

suppressPackageStartupMessages({
    library(DelayedArray)
    library(HDF5Array)
})

#' Absolute paths of every HDF5 file a DelayedArray reads from.
hdf5_seed_paths <- function(x) {
    if (!is(x, "DelayedArray")) return(character(0))
    p <- DelayedArray::seedApply(x, function(s) {
        if (is(s, "HDF5ArraySeed")) s@filepath else NA_character_
    })
    unique(stats::na.omit(unlist(p)))
}

#' Load a per-chromosome BSseq object, from whichever source is intact.
#'
#' @param region caudate | dlpfc | hippocampus
#' @param chrom chromosome, e.g. 22
#' @param require_smoothed stop unless the object carries usable smoothed fits
load_bsobj <- function(region, chrom, root = repo_root(),
                       require_smoothed = TRUE) {
    paths <- load_config("paths", root = root)
    tmpl <- paths$wgbs_bsobj_hdf5se[[region]]

    if (!is.null(tmpl)) {
        d <- gsub("{chrom}", as.character(chrom), tmpl, fixed = TRUE)
        if (!dir.exists(d)) {
            stop("HDF5SummarizedExperiment directory declared for ", region,
                 " but missing: ", d,
                 "\n  Set wgbs_bsobj_hdf5se.", region,
                 " to null in config/paths.yml to fall back to the .rda.")
        }
        message("[load] ", d, " (HDF5SummarizedExperiment)")
        bs <- HDF5Array::loadHDF5SummarizedExperiment(d)
    } else {
        f <- resolve_path("wgbs_bsobj_template", region = region, chrom = chrom,
                          root = root, check = TRUE)
        message("[load] ", f, " (.rda)")
        e <- new.env(parent = emptyenv())
        loaded <- load(f, envir = e)
        if (length(loaded) != 1) {
            stop("Expected exactly one object in ", f, ", found: ",
                 paste(loaded, collapse = ", "))
        }
        bs <- get(loaded, envir = e)
    }

    ## Fail here, in seconds, rather than deep inside an array job. Every assay
    ## the pipeline can touch must actually be readable -- the /tmp-backed `coef`
    ## that motivated this function passed every check except being read.
    need <- c("M", "Cov")
    if (require_smoothed) {
        if (!bsseq::hasBeenSmoothed(bs)) {
            stop("BSobj for ", region, " chr", chrom, " is not smoothed, but ",
                 "00_prepare.R calls getMeth(type = 'smooth').")
        }
        need <- c(need, "coef")
    }
    missing_assays <- setdiff(need, names(bs@assays@data))
    if (length(missing_assays) > 0) {
        stop("BSobj for ", region, " chr", chrom, " lacks assay(s): ",
             paste(missing_assays, collapse = ", "))
    }
    for (a in need) {
        probe <- try(as.matrix(bs@assays@data[[a]][seq_len(min(5L, nrow(bs))),
                                                   seq_len(min(2L, ncol(bs))),
                                                   drop = FALSE]),
                     silent = TRUE)
        if (inherits(probe, "try-error")) {
            stale <- hdf5_seed_paths(bs@assays@data[[a]])
            stop("Assay '", a, "' of ", region, " chr", chrom,
                 " is unreadable:\n  ",
                 conditionMessage(attr(probe, "condition")),
                 if (length(stale)) paste0("\n  HDF5 backing: ",
                                           paste(stale, collapse = ", ")) else "")
        }
    }
    bs
}
