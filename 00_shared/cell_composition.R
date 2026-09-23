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


## Where the two deconvolution modalities may be used.
##
## AGENTS.md 7.4 is asymmetric on purpose. RNA MuSiC adjustment is required in
## every region, unconditionally. DNAm scMD adjustment is required only "when the
## integration gate passes", and it passes in caudate alone: neuronal
## RNA-vs-DNAm concordance is rho 0.72 (FDR 2e-45) in caudate and is null in
## DLPFC (-0.035, p 0.66) and hippocampus (-0.103, p 0.27). An scMD-derived
## covariate outside caudate is therefore an adjustment for an estimate that
## does not track the composition it claims to measure.
##
## The thresholds are read from config/cell_deconvolution.yml:validation on every
## call rather than restated, and the gate is recomputed from the concordance
## table rather than trusted to a cached flag, so a re-deconvolution cannot leave
## a stale PASS behind. Moved here 2026-09-23 from
## 09b_aging_application/_h/age_functions.R, which defined it first; Module 04
## needs the identical decision and AGENTS.md 5.3 forbids a second copy.
scmd_gate_passes <- function(region, gate_table, root = repo_root()) {
    val <- load_config("cell_deconvolution", root = root)$validation
    tab <- data.table::fread(file.path(root, gate_table))
    ## Selection computed OUTSIDE `[`: inside it the bare name `region` would
    ## resolve to the column, matching every row (see gates.R).
    keep <- tab$region == region & tab$broad_class == "Total_neuron"
    row <- tab[which(keep)]
    if (nrow(row) != 1L) stop("No Total_neuron concordance row for ", region)
    isTRUE(row$rho >= as.numeric(val$min_neuronal_spearman_rho) &&
           row$neuron_fdr <= as.numeric(val$max_neuronal_fdr))
}


## Per-region proportion files and the scMD gate, resolved from config.
##
## `config/repeat_annotations.yml:covariate_sources.cell_composition_pcs` names
## `config/cell_deconvolution.yml` as the source of these estimates, so the
## templates and the gate table are taken from there instead of being spelled out
## in an analysis script (AGENTS.md 9: Quest paths live in configuration).
cell_composition_sources <- function(region, root = repo_root(),
                                     n_pcs = 3L) {
    cfg <- load_config("cell_deconvolution", root = root)
    fill <- function(tmpl) sub("{region}", region, tmpl, fixed = TRUE)
    gate_table <- file.path(cfg$paths$output_root,
                            "dnam-scmd-rna-concordance.tsv")
    music_f <- fill(cfg$paths$rna_proportions_template)
    scmd_f  <- fill(cfg$paths$dnam_proportions_template)
    scmd_ok <- scmd_gate_passes(region, gate_table, root = root)
    list(region = region,
         gate_table = gate_table,
         scmd_integration_gate = if (scmd_ok) "PASS" else "FAIL",
         music_file = music_f,
         scmd_file = scmd_f,
         n_pcs = as.integer(n_pcs),
         ## RNA MuSiC is unconditional; scMD is fitted only where the gate holds.
         pc_sets = c(list(music = cell_composition_pcs(
                              file.path(root, music_f), n_pcs = n_pcs)),
                     if (scmd_ok) list(scmd = cell_composition_pcs(
                              file.path(root, scmd_f), n_pcs = n_pcs))))
}


## Per-VMR cell-composition R^2, for one or more donor PC sets.
##
## The donor deconvolution estimates are per DONOR, not per VMR, so the per-VMR
## quantity is how strongly a VMR's methylation tracks composition: the R^2 of
## its methylation across donors on the donor cell-proportion PCs.
##
## Every PC set is evaluated against the SAME phenotype read, because a VMR's
## .phen file is one small read repeated ~11,000 times and reading it once per
## modality doubles the cost of the feature build for no gain.
vmr_composition_r2 <- function(phen_dir, chrom, start, end, pc_sets,
                               min_donors = 20L) {
    if (!length(pc_sets)) stop("vmr_composition_r2(): no PC sets given")
    n <- length(chrom)
    if (length(start) != n || length(end) != n) {
        stop("vmr_composition_r2(): chrom/start/end lengths differ")
    }
    out <- matrix(NA_real_, nrow = n, ncol = length(pc_sets),
                  dimnames = list(NULL, names(pc_sets)))
    for (i in seq_len(n)) {
        f <- file.path(phen_dir, sprintf("%s_%d_%d_meth.phen",
                                         chrom[i], start[i], end[i]))
        if (!file.exists(f)) next
        ph <- data.table::fread(f, header = FALSE,
                                col.names = c("FID", "IID", "y"),
                                colClasses = list(character = 1:2))
        for (s in names(pc_sets)) {
            pcs <- pc_sets[[s]]
            idx <- match(ph$FID, rownames(pcs))
            ok <- !is.na(idx) & is.finite(ph$y)
            if (sum(ok) < min_donors) next
            out[i, s] <- summary(stats::lm(
                ph$y[ok] ~ pcs[idx[ok], , drop = FALSE]))$r.squared
        }
    }
    out
}
