#### Chromosome and chunk ordering (v2 revision) ####
##
## Defects addressed here:
##
## V4  Legacy per-chromosome arrays ran 1-24 over c(1:22, "X", "Y") and merged
##     the results into one catalog. The C->T SNP masks only exist for autosomes,
##     so chrX/chrY CpGs were never masked -- caudate ended up with a 3x excess
##     of sex-chromosome VMRs (431 vs 143/147). v2 keeps sex chromosomes out of
##     the primary catalog and writes them to an explicit excluded/ manifest.
##
## V7  res_cpg_meth.phen columns were assembled in list.files() order, which is
##     lexicographic: chunk_10000 sorts before chunk_5000. The stored column
##     order therefore did not match genomic order. sort_chunks() sorts on the
##     numeric index embedded in the filename instead.
##
## Also V7: the legacy combine step cbind()ed chunks while assuming identical
## FID/IID order across them, without checking. verify_fid_iid() checks.

suppressPackageStartupMessages({
    library(data.table)
})

#' Chromosomes for a run, in deterministic numeric order.
#'
#' @param include_sex override the config policy (used by the X/Y export step)
#' @return character vector, autosomes ascending then X, Y if requested
chrom_order <- function(include_sex = NULL, root = repo_root()) {
    th <- load_config("thresholds", root = root)
    auto <- as.character(sort(as.integer(th$chromosomes$primary)))
    if (is.null(include_sex)) include_sex <- isTRUE(th$chromosomes$include_sex_in_primary)
    if (!include_sex) return(auto)
    c(auto, sex_chroms(th))
}

#' Sex chromosome names from config, guarding the YAML 1.1 boolean trap.
#'
#' A bare `Y` in YAML 1.1 parses as logical TRUE, which silently turns
#' `sex: [X, Y]` into c("X", "TRUE") and drops chrY from every sex-chromosome
#' check. The config quotes them; this catches it if someone unquotes them again.
sex_chroms <- function(th = load_config("thresholds")) {
    raw <- th$chromosomes$sex
    if (any(vapply(raw, is.logical, logical(1)))) {
        stop("config/thresholds.yml chromosomes.sex contains a boolean. ",
             "YAML 1.1 reads a bare Y as true -- quote the values: sex: [\"X\", \"Y\"]")
    }
    as.character(unlist(raw))
}

#' Is this chromosome part of the primary catalog?
is_primary_chrom <- function(chrom, root = repo_root()) {
    as.character(chrom) %in% chrom_order(include_sex = FALSE, root = root)
}

#' Does a C->T SNP mask exist for this chromosome?
#'
#' Sex chromosomes have no mask file. The legacy code handled this with
#' `if (!chr %in% c("X","Y"))` and then merged the unmasked results anyway.
has_ct_mask <- function(chrom, root = repo_root()) {
    th <- load_config("thresholds", root = root)
    as.character(chrom) %in% as.character(th$chromosomes$ct_mask_available)
}

#' Sort chunk files by their embedded numeric start index, not lexicographically.
#'
#' Filenames look like cpg_meth_1_5000.tsv, cpg_meth_5001_10000.tsv, ...
#' list.files() would order 10001_15000 before 5001_10000.
sort_chunks <- function(files) {
    if (length(files) == 0) return(files)
    idx <- vapply(basename(files), function(f) {
        nums <- as.numeric(regmatches(f, gregexpr("[0-9]+", f))[[1]])
        if (length(nums) == 0) NA_real_ else nums[1]
    }, numeric(1))
    if (anyNA(idx)) {
        stop("Could not extract a numeric index from chunk filename(s): ",
             paste(basename(files)[is.na(idx)], collapse = ", "))
    }
    assert_no_dups(idx, "chunk start indices")
    files[order(idx)]
}

#' Verify every chunk carries identical FID/IID in identical order before cbind.
#'
#' @param chunks list of data.tables, each with leading FID and IID columns
verify_fid_iid <- function(chunks, labels = NULL) {
    if (length(chunks) == 0) stop("verify_fid_iid(): no chunks given")
    if (is.null(labels)) labels <- paste0("chunk_", seq_along(chunks))

    for (i in seq_along(chunks)) {
        if (!all(c("FID", "IID") %in% names(chunks[[i]]))) {
            stop("Missing FID/IID columns in ", labels[[i]])
        }
    }

    ref <- chunks[[1]][, .(FID = as.character(FID), IID = as.character(IID))]
    for (i in seq_along(chunks)[-1]) {
        this <- chunks[[i]][, .(FID = as.character(FID), IID = as.character(IID))]
        if (!identical(ref$FID, this$FID) || !identical(ref$IID, this$IID)) {
            n_diff <- if (nrow(ref) == nrow(this)) sum(ref$FID != this$FID) else NA
            stop("Donor rows differ between ", labels[[1]], " and ", labels[[i]],
                 ": ", nrow(ref), " vs ", nrow(this), " rows",
                 if (!is.na(n_diff)) paste0(", ", n_diff, " mismatched FIDs"),
                 ".\n  Chunks may not be cbind()ed unless donor order is identical.")
        }
    }
    invisible(ref)
}

#' Sort a VMR/CpG table into genomic order (numeric chromosome, then position).
sort_genomic <- function(dt, chr_col = "chr", pos_col = "start") {
    dt <- data.table::as.data.table(dt)
    chr_clean <- sub("^chr", "", as.character(dt[[chr_col]]))
    ## X -> 23, Y -> 24 so ordering is total and deterministic
    chr_num <- suppressWarnings(as.integer(chr_clean))
    chr_num[chr_clean == "X"] <- 23L
    chr_num[chr_clean == "Y"] <- 24L
    if (anyNA(chr_num)) {
        stop("Unrecognized chromosome value(s): ",
             paste(unique(chr_clean[is.na(chr_num)]), collapse = ", "))
    }
    dt[order(chr_num, dt[[pos_col]])]
}
