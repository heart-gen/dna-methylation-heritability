#### 09b_aging_application -- module-local functions ####
##
## Sourced by 01, 02 and 05. Everything here is linear algebra on a donor x VMR
## methylation matrix, because the donor bootstraps refit every VMR a thousand
## times: one QR per design, applied to all VMRs at once,
## replaces thousands of lm() calls with a few matrix products. The smoke test
## checks that these agree with per-VMR lm() to 1e-10.

## Read one Module 01 phenotype per eligible VMR into a donor x VMR matrix.
##
## Covariates come from 00_shared/locus_io.R::load_locus_phenotype(), the reader
## Module 02 used. Module 01 writes one covariate file per chromosome; they are
## asserted identical to the first chromosome's, so a single donor table is
## correct for every VMR. Each phenotype is ALIGNED to the reference donor order
## by ID, never assumed to share it (AGENTS.md 7.1).
load_age_inputs <- function(tasks, vmr_run_dir, catalog_cohort, covar_prefix) {
    ref <- NULL
    cov_seen <- list()
    Y <- matrix(NA_real_, nrow = 0L, ncol = nrow(tasks))
    for (i in seq_len(nrow(tasks))) {
        task <- as.list(tasks[i])
        ph <- load_locus_phenotype(task = task, vmr_run_dir = vmr_run_dir,
                                   cohort = catalog_cohort,
                                   covar_prefix = covar_prefix)
        if (!is.null(ph$pcs)) {
            stop("Estimation-cell genotype PCs present for ", task$vmr_id,
                 "; this module is prespecified for the AA arm only.")
        }
        chrom <- as.character(task$chrom)
        if (is.null(ref)) {
            ref <- merge(ph$covar, ph$qcovar, by = c("FID", "IID"), all = FALSE)
            ref <- ref[order(ref$FID), , drop = FALSE]
            if (anyDuplicated(ref$FID)) stop("Duplicate donors in Module 01 covariates")
            if (nrow(ref) != nrow(ph$covar) || nrow(ref) != nrow(ph$qcovar)) {
                stop("covar and qcovar donor sets differ for ", chrom)
            }
            Y <- matrix(NA_real_, nrow = nrow(ref), ncol = nrow(tasks),
                        dimnames = list(ref$FID, tasks$vmr_id))
        }
        if (is.null(cov_seen[[chrom]])) {
            this <- merge(ph$covar, ph$qcovar, by = c("FID", "IID"), all = FALSE)
            this <- this[order(this$FID), , drop = FALSE]
            rownames(this) <- NULL
            chk <- ref; rownames(chk) <- NULL
            if (!isTRUE(all.equal(this, chk, check.attributes = FALSE))) {
                stop("Module 01 covariates differ between chromosomes (", chrom,
                     "); one donor table cannot serve every VMR.")
            }
            cov_seen[[chrom]] <- TRUE
        }
        if (anyDuplicated(ph$phenotype$FID)) {
            stop("Duplicate donors in phenotype for ", task$vmr_id)
        }
        idx <- match(ref$FID, ph$phenotype$FID)
        if (anyNA(idx) || nrow(ph$phenotype) != nrow(ref)) {
            stop("Phenotype donors for ", task$vmr_id, " differ from the ",
                 "covariate donors (AGENTS.md 10.1: missing donors fail loudly).")
        }
        Y[, i] <- as.numeric(ph$phenotype$phenotype[idx])
    }
    list(Y = Y, donors = data.table::as.data.table(ref))
}

## Donor-level design for one spec. Returns the rows (donor indices into the
## full table), the model matrix, and the permutation strata. Age is in decades.
build_spec_design <- function(spec, donors, cfg, cellpcs = list()) {
    ac <- cfg$age_model
    unit <- as.numeric(ac$age_unit_years)
    ctrl <- as.character(ac$diagnosis_control_label)
    keep <- rep(TRUE, nrow(donors))
    if (identical(spec$donors, "controls")) keep <- keep & donors$diagnosis == ctrl
    if (!is.null(spec$min_age)) keep <- keep & donors$age >= as.numeric(spec$min_age)

    d <- data.frame(age = donors$age / unit,
                    sex = factor(donors$sex),
                    diagnosis = factor(donors$diagnosis))
    extra <- as.character(unlist(spec$extra_terms))
    for (term in extra) {
        src <- sub("_cellPC[0-9]+$", "", term)
        pcs <- cellpcs[[src]]
        if (is.null(pcs)) stop("Spec needs ", src, " cell PCs, which were not built")
        col <- sub("^.*_", "", term)
        col <- sub("cellPC", "PC", col)
        if (!col %in% colnames(pcs)) stop("No ", col, " in ", src, " cell PCs")
        v <- pcs[match(donors$FID, rownames(pcs)), col]
        d[[term]] <- v
        keep <- keep & is.finite(v)
    }
    terms <- setdiff(c("age", "sex", "diagnosis"), unlist(spec$drop_terms))
    terms <- c(terms, extra)
    d <- d[keep, , drop = FALSE]
    for (f in c("sex", "diagnosis")) if (f %in% terms) d[[f]] <- droplevels(d[[f]])
    X <- stats::model.matrix(stats::reformulate(terms), data = d)
    if (qr(X)$rank < ncol(X)) {
        stop("Design for spec is rank-deficient: ", paste(colnames(X), collapse = ","))
    }
    list(rows = which(keep), X = X, age_col = which(colnames(X) == "age"),
         strata = as.character(donors$diagnosis[keep]),
         n = sum(keep),
         n_controls = sum(donors$diagnosis[keep] == ctrl),
         n_cases = sum(donors$diagnosis[keep] != ctrl))
}

## Age coefficient, SE, t and p for every VMR at once. Y is n x V.
## The generic matrix fitter lives in 00_shared/axis_inference.R, which Module 10
## also uses; this keeps 09b's call sites and its age-specific naming.
fit_age_matrix <- function(X, Y, age_col) fit_scalar_matrix(X, Y, age_col)

## resample_rows(), prepare_axis(), axis_estimate(), block_jackknife_se()
## and combined_inference() moved to 00_shared/axis_inference.R on
## 2026-09-19, arithmetic unchanged, when Module 10 needed the same
## machinery (AGENTS.md 5.3: refactor shared logic rather than duplicate).
## 09b/tests/ checks the shared versions against the definitions this file
## carried at commit abb0c6789.

## Per-VMR outcomes from one set of age fits.
##
## debiased_sq_beta = beta_hat^2 - SE^2 is UNBIASED for beta^2 whatever the SE.
## That is the property the design rests on: a VMR with strong local SNP control
## carries its genetic variance in the age model's residual, so its SE is
## larger, and any outcome whose expectation depends on SE (|beta_hat|, its
## rank, |t|, -log10 p) couples to the score mechanically. See config/aging.yml.
age_outcomes <- function(fit) {
    list(debiased_sq_beta = fit$beta^2 - fit$se^2,
         signed_beta_age = fit$beta)
}


## scMD integration gate, re-read from the concordance table on every run
## (config/cell_deconvolution.yml:validation).
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
