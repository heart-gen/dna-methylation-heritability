#### Donor identity and row alignment (v2 revision) ####
##
## This file exists because of defect V1.
##
## The legacy code (vmr-analysis/*/_h/02b.res_var.R:39-40) did:
##
##     meth_levels <- meth_levels[match(valid_ids, brain_id), , drop = FALSE]
##     pc_filt     <- pc[V1 %in% valid_ids, -1, with = FALSE]
##
## The response matrix is REORDERED to `valid_ids` order. The design matrix is
## merely SUBSET, so it keeps pc.csv's original order. Unless those two orders
## happen to coincide, every donor's methylation is regressed against a
## different donor's principal components. Measured impact: caudate chr1 had
## 153/153 rows wrong, DLPFC 92/96, hippocampus 98/101. The 99th-percentile SD
## cutoff shifted +9.2% and ~10% of seed CpGs changed (Jaccard 0.82), which
## means every VMR set in the repository is invalid.
##
## The bug is subtle because `%in%` subsetting looks like alignment and produces
## a matrix of exactly the right dimensions. It fails silently and forever.
##
## The rule in v2: never subset two things and assume they line up. Both sides
## go through align_by_id(), which reorders BOTH by match() and then asserts the
## resulting ID vectors are identical.

suppressPackageStartupMessages({
    library(data.table)
})

#' Align two tables to a common, explicitly ordered set of donor IDs.
#'
#' Reorders both inputs by `match()` against `ids` -- never by `%in%` -- and
#' verifies the result before returning. Any ID in `ids` that is missing from
#' either side is an error, not a silently dropped row.
#'
#' @param x,y data.frame/data.table/matrix, one row per donor
#' @param id_x,id_y donor ID vectors, parallel to rows of x and y
#' @param ids the ordered ID set to align to; defaults to the intersection
#'   taken in id_x order
#' @return list(x, y, ids) with nrow(x) == nrow(y) == length(ids)
align_by_id <- function(x, y, id_x, id_y, ids = NULL) {
    id_x <- as.character(id_x)
    id_y <- as.character(id_y)

    if (nrow_any(x) != length(id_x)) {
        stop("align_by_id(): x has ", nrow_any(x), " rows but id_x has ",
             length(id_x), " entries")
    }
    if (nrow_any(y) != length(id_y)) {
        stop("align_by_id(): y has ", nrow_any(y), " rows but id_y has ",
             length(id_y), " entries")
    }

    assert_no_dups(id_x, "id_x")
    assert_no_dups(id_y, "id_y")

    if (is.null(ids)) {
        ids <- id_x[id_x %in% id_y]
    } else {
        ids <- as.character(ids)
        assert_no_dups(ids, "ids")
        missing_x <- setdiff(ids, id_x)
        missing_y <- setdiff(ids, id_y)
        if (length(missing_x) > 0 || length(missing_y) > 0) {
            stop("align_by_id(): requested IDs absent from input.",
                 if (length(missing_x)) paste0("\n  missing from x: ",
                     paste(head(missing_x, 10), collapse = ", ")),
                 if (length(missing_y)) paste0("\n  missing from y: ",
                     paste(head(missing_y, 10), collapse = ", ")))
        }
    }

    if (length(ids) == 0) {
        stop("align_by_id(): no donors in common between x and y")
    }

    ix <- match(ids, id_x)
    iy <- match(ids, id_y)
    if (anyNA(ix) || anyNA(iy)) {
        stop("align_by_id(): internal match failure -- this should be ",
             "unreachable and indicates duplicate or NA IDs")
    }

    x_aligned <- subset_rows(x, ix)
    y_aligned <- subset_rows(y, iy)

    ## The assertion the legacy code lacked.
    stopifnot(identical(as.character(id_x[ix]), ids),
              identical(as.character(id_y[iy]), ids),
              nrow_any(x_aligned) == length(ids),
              nrow_any(y_aligned) == length(ids))

    list(x = x_aligned, y = y_aligned, ids = ids)
}

nrow_any <- function(z) if (is.null(dim(z))) length(z) else nrow(z)

subset_rows <- function(z, i) {
    if (is.null(dim(z))) return(z[i])
    if (data.table::is.data.table(z)) return(z[i, , drop = FALSE])
    z[i, , drop = FALSE]
}

#' Read a PLINK .psam, with or without a leading '#FID IID SEX' header.
#'
#' Defect V8: the AA .psam is headerless while the all_individuals .psam is not.
#' Reading the AA file with header = TRUE consumed its first data row and
#' silently dropped donor Br2585 from every downstream analysis. Reading
#' headerless and discarding any '#'-prefixed row handles both files.
read_psam <- function(psam_file) {
    if (!file.exists(psam_file)) stop("psam not found: ", psam_file)
    samples <- data.table::fread(psam_file, header = FALSE,
                                 colClasses = "character")
    if (ncol(samples) < 2) {
        stop("psam has fewer than 2 columns: ", psam_file)
    }
    samples <- samples[, 1:2]
    data.table::setnames(samples, c("FID", "IID"))
    samples <- samples[!startsWith(FID, "#")]
    if (nrow(samples) == 0) stop("psam contained no donor rows: ", psam_file)
    if (anyNA(samples$FID) || anyNA(samples$IID)) {
        stop("NA IDs in psam: ", psam_file)
    }
    assert_no_dups(samples$FID, paste0("FID in ", basename(psam_file)))
    samples[]
}

#' Fail loudly on duplicate identifiers (AGENTS.md 10.1).
assert_no_dups <- function(ids, what = "identifiers") {
    dups <- unique(ids[duplicated(ids)])
    if (length(dups) > 0) {
        stop("Duplicate ", what, ": ", paste(head(dups, 10), collapse = ", "),
             if (length(dups) > 10) paste0(" (and ", length(dups) - 10, " more)"))
    }
    invisible(TRUE)
}

#' Fail loudly on missing identifiers.
assert_present <- function(ids, required, what = "donors") {
    missing <- setdiff(as.character(required), as.character(ids))
    if (length(missing) > 0) {
        stop("Missing ", what, ": ", paste(head(missing, 10), collapse = ", "),
             if (length(missing) > 10) paste0(" (and ", length(missing) - 10, " more)"))
    }
    invisible(TRUE)
}

#' Compare the observed analysis-set size against the PI-locked design count.
#'
#' AGENTS.md 14: "sample counts differ from the design without explanation" is a
#' stop condition. While design_n is null (not yet locked) we record and warn
#' rather than stop, so the first run can establish the number.
assert_expected_n <- function(observed_n, cohort, region, root = repo_root()) {
    cohorts <- load_config("cohorts", root = root)
    entry <- cohorts$donor_counts[[cohort]][[region]]
    if (is.null(entry)) {
        stop("No donor_counts entry for ", cohort, "/", region,
             " in config/cohorts.yml")
    }

    if (!is.null(entry$phenotype_n) && observed_n > entry$phenotype_n) {
        stop("Observed n (", observed_n, ") exceeds the phenotype-table ceiling (",
             entry$phenotype_n, ") for ", cohort, "/", region,
             ". The analysis set cannot be larger than the donors the phenotype ",
             "table supplies -- check the region and race filters.")
    }

    if (is.null(entry$design_n)) {
        message("[design] ", cohort, "/", region, ": observed n = ", observed_n,
                " (design_n not yet locked; ceiling ", entry$phenotype_n, "). ",
                "Record this in config/cohorts.yml and lock it with the PI ",
                "before production.")
        return(invisible(observed_n))
    }

    if (observed_n != entry$design_n) {
        stop("Sample count mismatch for ", cohort, "/", region, ": observed ",
             observed_n, ", design ", entry$design_n,
             ".\n  AGENTS.md 14: stop and request PI direction rather than ",
             "proceeding with an unexplained count.")
    }
    invisible(observed_n)
}

#' Checksum of the ORDERED donor list (AGENTS.md 9).
#'
#' Order matters: two runs over the same donors in different orders are not the
#' same run, and a changed checksum is how a downstream module detects that its
#' upstream inputs moved.
donor_checksum <- function(ids) {
    s <- paste(as.character(ids), collapse = "\n")
    if (requireNamespace("digest", quietly = TRUE)) {
        return(digest::digest(s, algo = "sha256", serialize = FALSE))
    }
    tf <- tempfile(); on.exit(unlink(tf), add = TRUE)
    writeLines(s, tf)
    out <- tryCatch(system2("sha256sum", shQuote(tf), stdout = TRUE),
                    error = function(e) NA_character_)
    if (length(out) == 0 || is.na(out[1])) return(NA_character_)
    sub(" .*$", "", out[1])
}
