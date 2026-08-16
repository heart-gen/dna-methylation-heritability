#### WGBS BSseq object loading (v2 revision) ####
##
## The per-chromosome BSobj .rda files are thin DelayedArray views: the actual
## methylation counts live in a genome-wide HDF5 file, referenced by an absolute
## path baked into each object when it was saved. Those baked-in paths have gone
## stale (V14), so loading a BSobj is not a plain load() and every v2 module goes
## through load_bsobj().
##
## AGENTS.md 9: Quest paths belong in configuration. The authoritative HDF5
## location per region is config/paths.yml:wgbs_hdf5_assays, and this file
## repoints the object at it -- loudly, never silently.

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

#' Repoint every HDF5 seed of a DelayedArray at `new_path`.
repoint_hdf5 <- function(x, new_path) {
    DelayedArray::modify_seeds(x, function(s) {
        if (is(s, "HDF5ArraySeed")) s@filepath <- new_path
        s
    })
}

#' Load a per-chromosome BSseq object with its HDF5 backing repaired.
#'
#' V14: the DLPFC BSobj .rda files reference
#'   .../new-data/dlpfc/_m/combined_hdf5/dlpfc_assays.h5
#' but the file on disk is `assays.h5` -- the data owner renamed it. Every
#' 00_prepare.R task for DLPFC died on the missing file. Verified 2026-08-15
#' that the two are the same data: identical dims (29,401,795 x 176), identical
#' donor names in identical order, and Cov values bit-identical to the
#' genome-wide se.rds for the chr22 rows checked. Hippocampus and caudate
#' already resolve correctly and are left untouched by the repair.
#'
#' The `coef` assay (smoothing output) points at a /tmp scratch directory from
#' the session that produced it and is unrecoverable for every region. v2 never
#' reads it, so it is dropped rather than left as a landmine that fails deep
#' inside an array job.
#'
#' @param region caudate | dlpfc | hippocampus
#' @param chrom chromosome, e.g. 22
#' @param assays_used assays v2 requires; the rest are dropped
load_bsobj <- function(region, chrom, root = repo_root(),
                       assays_used = c("M", "Cov")) {
    f <- resolve_path("wgbs_bsobj_template", region = region, chrom = chrom,
                      root = root, check = TRUE)
    message("[load] ", f)
    e <- new.env(parent = emptyenv())
    loaded <- load(f, envir = e)
    if (length(loaded) != 1) {
        stop("Expected exactly one object in ", f, ", found: ",
             paste(loaded, collapse = ", "))
    }
    bs <- get(loaded, envir = e)

    missing_assays <- setdiff(assays_used, names(bs@assays@data))
    if (length(missing_assays) > 0) {
        stop("BSobj for ", region, " chr", chrom, " lacks assay(s): ",
             paste(missing_assays, collapse = ", "))
    }
    ## Drop everything we do not use, notably the /tmp-backed `coef`.
    bs@assays@data <- bs@assays@data[assays_used]

    h5 <- resolve_path(paste0("wgbs_hdf5_assays.", region), root = root,
                       check = TRUE)
    for (a in assays_used) {
        cur <- hdf5_seed_paths(bs@assays@data[[a]])
        if (length(cur) == 0) next            # not HDF5-backed; nothing to repair
        if (length(cur) > 1) {
            stop("Assay '", a, "' of ", region, " chr", chrom,
                 " spans multiple HDF5 files, which v2 cannot repoint safely:\n  ",
                 paste(cur, collapse = "\n  "))
        }
        if (cur == h5) next                   # already correct (caudate, hippocampus)
        if (file.exists(cur)) {
            ## Two live files disagreeing is a provenance question, not ours to
            ## resolve by preferring one.
            stop("Assay '", a, "' of ", region, " chr", chrom, " points at\n  ",
                 cur, "\nwhich EXISTS but is not the configured backing store\n  ",
                 h5, "\nRefusing to guess which is authoritative.")
        }
        message("[hdf5] ", region, " chr", chrom, " assay '", a,
                "': stale reference ", basename(cur), " -> ", h5, " (V14)")
        bs@assays@data[[a]] <- repoint_hdf5(bs@assays@data[[a]], h5)
    }

    ## Prove the repair before an array job spends hours on it.
    for (a in assays_used) {
        probe <- try(as.matrix(bs@assays@data[[a]][seq_len(min(5L, nrow(bs))),
                                                   seq_len(min(2L, ncol(bs))),
                                                   drop = FALSE]),
                     silent = TRUE)
        if (inherits(probe, "try-error")) {
            stop("Assay '", a, "' of ", region, " chr", chrom,
                 " is unreadable after HDF5 repair:\n  ", conditionMessage(attr(probe, "condition")))
        }
    }
    bs
}
