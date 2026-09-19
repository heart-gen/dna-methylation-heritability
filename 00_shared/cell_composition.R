## Donor cell-composition principal components.
##
## Lifted from 04_repeat_repressive_architecture/_h/01_build_features.R, where
## it built the donor PCs behind `cell_composition_r2`, so that Module 04 and
## 09b_aging_application construct composition covariates identically rather
## than from two copies that could drift (AGENTS.md 5.3).
##
## `prop_file` is a long table (sample_id, cell_type, proportion), as written by
## inputs/cell_proportions for both RNA MuSiC and DNAm scMD. Cell types with
## zero or undefined variance are dropped before the PCA, and the PC count is
## capped at ncol - 1 because proportions close to one: the last PC of a
## compositional matrix carries only rounding.
##
## The PCA runs over EVERY donor in the file, not only the caller's analysis
## donors. That is what Module 04 did and is kept so its accepted numbers
## reproduce; callers subset rows afterwards by rowname (the donor brnum).
cell_composition_pcs <- function(prop_file, n_pcs = 3L) {
    if (!file.exists(prop_file)) {
        stop("Cell-proportion estimates not found: ", prop_file)
    }
    props <- data.table::fread(prop_file)
    need <- c("sample_id", "cell_type", "proportion")
    if (!all(need %in% names(props))) {
        stop("Cell-proportion table lacks ", paste(setdiff(need, names(props)),
             collapse = ", "), ": ", prop_file)
    }
    wide <- data.table::dcast(props, sample_id ~ cell_type,
                              value.var = "proportion")
    pmat <- as.matrix(wide[, -1])
    rownames(pmat) <- wide$sample_id
    pmat <- pmat[, apply(pmat, 2, function(z) is.finite(stats::var(z)) &&
                                              stats::var(z) > 0), drop = FALSE]
    n_pc <- min(as.integer(n_pcs), ncol(pmat) - 1L)
    if (n_pc < 1L) stop("No usable cell-type columns in ", prop_file)
    stats::prcomp(pmat, center = TRUE, scale. = TRUE)$x[, seq_len(n_pc),
                                                        drop = FALSE]
}
