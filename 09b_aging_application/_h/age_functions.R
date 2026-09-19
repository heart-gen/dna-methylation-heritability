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
fit_age_matrix <- function(X, Y, age_col) {
    q <- qr(X)
    beta_all <- qr.coef(q, Y)
    resid <- qr.resid(q, Y)
    df <- nrow(X) - q$rank
    sigma2 <- colSums(resid^2) / df
    xtx_inv <- chol2inv(qr.R(q))
    ## qr() may pivot; map the age column through the pivot.
    piv <- q$pivot
    a <- which(piv == age_col)
    se <- sqrt(sigma2 * xtx_inv[a, a])
    beta <- if (is.matrix(beta_all)) beta_all[age_col, ] else beta_all[age_col]
    t <- beta / se
    list(beta = as.numeric(beta), se = as.numeric(se), t = as.numeric(t),
         p = as.numeric(2 * stats::pt(-abs(t), df)), df = df)
}

## Resample donors with replacement WITHIN strata (diagnosis), so every draw
## keeps the observed case/control counts and the design stays estimable.
resample_rows <- function(strata) {
    idx <- integer(0)
    for (s in unique(strata)) {
        pool <- which(strata == s)
        idx <- c(idx, pool[sample.int(length(pool), length(pool), replace = TRUE)])
    }
    idx
}

## Axis design for a fixed set of VMR rows. Only the outcome changes between
## bootstrap draws, so the QR is taken once and every refit is a single qr.coef.
prepare_axis <- function(dt, predictor, covariates, rows = NULL) {
    cols <- c(predictor, covariates)
    ok <- stats::complete.cases(dt[, cols, with = FALSE])
    if (!is.null(rows)) ok <- ok & rows
    Z <- cbind(`(Intercept)` = 1, as.matrix(dt[ok, cols, with = FALSE]))
    storage.mode(Z) <- "double"
    q <- qr(Z)
    if (q$rank < ncol(Z)) {
        stop("Axis design is rank-deficient for covariates: ",
             paste(covariates, collapse = ","))
    }
    list(ok = ok, Z = Z, qr = q, j = which(colnames(Z) == predictor),
         predictor = predictor, covariates = covariates)
}

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

## The axis estimate for one outcome on a subset of the axis rows. A
## "relative_to_mean" outcome is divided by its mean over the SAME rows, so the
## coefficient reads as the proportional change in mean squared age effect per
## SD of score -- unit-free, and so comparable between regions whose absolute
## age-effect scales differ. `drop` removes rows (the jackknife).
axis_estimate <- function(ax, y, scale = c("raw", "relative_to_mean"),
                          drop = NULL) {
    scale <- match.arg(scale)
    yy <- y[ax$ok]
    Z <- ax$Z
    if (!is.null(drop)) {
        keep <- !drop[ax$ok]
        yy <- yy[keep]; Z <- Z[keep, , drop = FALSE]
        q <- qr(Z)
    } else {
        q <- ax$qr
    }
    if (scale == "relative_to_mean") yy <- yy / mean(yy)
    unname(qr.coef(q, yy)[ax$j])
}

## Delete-one-chromosome weighted block jackknife (Busing, Meijer & van der
## Leeden 1999) -- the construction Modules 06 and 08 use. Weighted because
## chromosomes differ several-fold in VMR count. This is the VMR-level half of
## the variance: the VMRs are a sample of loci, and nearby loci are not
## independent, so the block is the chromosome.
block_jackknife_se <- function(ax, y, chrom, scale) {
    blocks <- sort(unique(chrom[ax$ok]))
    if (length(blocks) < 2L) return(NA_real_)
    full <- axis_estimate(ax, y, scale)
    n_tot <- sum(ax$ok)
    hj <- n_tot / vapply(blocks, function(b) sum(chrom[ax$ok] == b), numeric(1))
    theta <- vapply(blocks, function(b) axis_estimate(ax, y, scale,
                                                      drop = chrom == b), numeric(1))
    pseudo <- hj * full - (hj - 1) * theta
    sqrt(sum((pseudo - mean(pseudo))^2 / (hj - 1)) / length(blocks))
}

## Combined inference. Donor-bootstrap variance carries the donor-level half
## (shared donors correlate every VMR's estimate); the block jackknife carries
## the VMR-level half. Adding them double-counts the part of the per-VMR noise
## both see, which errs conservative; in simulation (tests/) the donor
## bootstrap alone rejected 57% of true nulls, the sum 0-3%.
combined_inference <- function(estimate, boot, se_jk, alpha = 0.05) {
    boot <- boot[is.finite(boot)]
    se <- sqrt(stats::var(boot) + se_jk^2)
    z <- stats::qnorm(1 - alpha / 2)
    list(se = se, se_bootstrap = stats::sd(boot), se_jackknife = se_jk,
         z = estimate / se, p = 2 * stats::pnorm(-abs(estimate / se)),
         ci_lower = estimate - z * se, ci_upper = estimate + z * se,
         n_bootstrap_used = length(boot))
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
