#### 02c_cis_greml_sensitivity -- module functions ####
##
## The GCTA calls themselves are 00_shared/gcta.R, shared with 02b. This file
## holds only what is particular to 02c: config, interpretation guards, unit
## IDs, and the chromosome-block jackknife for the aggregate statistics.

suppressPackageStartupMessages(library(data.table))

MODULE <- "02c_cis_greml_sensitivity"

load_cgs_config <- function(run_dir = NULL) {
    if (!is.null(run_dir)) {
        f <- file.path(run_dir, "code", "config", "cis_greml_sensitivity.yml")
        if (file.exists(f)) return(yaml::read_yaml(f))
    }
    load_config("cis_greml_sensitivity")
}

## The PI's three boundaries (2026-10-03), refused rather than defaulted.
assert_cgs_interpretation <- function(cfg) {
    for (nm in c("replaces_module_02_score", "per_vmr_absolute_h2_reportable",
                 "greml_significance_class_allowed",
                 "absolute_pve_interpretation_allowed",
                 "cross_region_level_comparison_allowed")) {
        if (!identical(cfg$interpretation[[nm]], FALSE)) {
            stop("config/cis_greml_sensitivity.yml interpretation.", nm,
                 " must be false: 02c validates the Module 02 ordering, it does ",
                 "not replace it, report per-VMR h2, or classify VMRs.")
        }
    }
    invisible(TRUE)
}

run_dir_for <- function(run_id) file.path(repo_root(), MODULE, "_m", "runs", run_id)

read_run_manifest <- function(run_dir) {
    m <- fread(file.path(run_dir, "manifest.tsv"), colClasses = "character")
    stats::setNames(m$value, m$field)
}

write_status <- function(rows, path) {
    rows <- as.data.table(rows)
    stopifnot(all(c("unit", "status", "reason") %in% names(rows)))
    write_atomic(rows, path)
}

unit_id <- function(vmr_id, mode) paste0(vmr_id, "|", mode)

## Flags carried on every emitted row.
cgs_flags <- function(dt) {
    dt[, `:=`(replaces_module_02_score = FALSE,
              per_vmr_absolute_h2_reportable = FALSE,
              greml_significance_class_allowed = FALSE,
              absolute_pve_interpretation_allowed = FALSE)]
    dt
}

## Delete-one-chromosome weighted block jackknife (Busing, Meijer & van der
## Leeden 1999) for an arbitrary statistic of a table -- the construction
## Modules 06, 08, 09b and 10 use. Loci on one chromosome are not independent,
## so the chromosome is the block. Returns the full-data estimate, the
## jackknife SE and a normal 95% interval.
block_jackknife_stat <- function(dt, block_col, stat_fn) {
    full <- stat_fn(dt)
    blocks <- sort(unique(dt[[block_col]]))
    if (length(blocks) < 2L || !is.finite(full)) {
        return(data.table(estimate = full, se = NA_real_, ci_low = NA_real_,
                          ci_high = NA_real_, n_blocks = length(blocks)))
    }
    n_tot <- nrow(dt)
    hj <- n_tot / vapply(blocks, function(b) sum(dt[[block_col]] == b), numeric(1))
    ## The row mask is built OUTSIDE data.table's `[`: inside it, a loop
    ## variable that shares a name with a column resolves to the column, and
    ## every deletion then silently drops all rows.
    blk <- dt[[block_col]]
    theta <- vapply(blocks, function(drop_block) {
        keep <- blk != drop_block
        stat_fn(dt[keep])
    }, numeric(1))
    pseudo <- hj * full - (hj - 1) * theta
    se <- sqrt(sum((pseudo - mean(pseudo))^2 / (hj - 1)) / length(blocks))
    data.table(estimate = full, se = se, ci_low = full - 1.96 * se,
               ci_high = full + 1.96 * se, n_blocks = length(blocks))
}
