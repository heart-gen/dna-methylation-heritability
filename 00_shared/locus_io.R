## Shared observed-locus reader.
##
## Stage 01 and the observed-regime simulation grid must consume byte-identical
## genotype, covariate and donor-alignment logic, otherwise the grid would
## characterise a different feature distribution than the one production
## actually produced. This file is the single implementation; neither caller
## reimplements any part of it.

## Returns a list with `status` one of "ok", "qc_failed" or "excluded". On "ok"
## the list carries genotype (donor-aligned, QC-filtered), phenotype y,
## covariates, and the source paths. Callers that simulate a phenotype ignore
## `y` but keep everything else, so n, num_snps, LD and p_eff are the observed
## values by construction.
## `apply_snp_qc = FALSE` returns the in-window genotype matrix WITHOUT the
## MAF/missingness filter, for callers that must perform every data-dependent
## genotype filter inside an outer training fold (AGENTS.md 7.3 names MAF,
## missingness and zero-variance filters among the things held-out donors must
## not influence). The eligibility gate still uses the QC'd variant count, so
## the task universe is identical either way -- only the returned matrix
## differs. Module 02 uses the default and is unaffected.
load_observed_locus <- function(task, cohort, vmr_run_dir, min_cis_variants,
                                expected_n = NA_integer_,
                                backing_tag = "lgv",
                                apply_snp_qc = TRUE) {
    chromosome_label <- sub("^chr", "", task$chrom, ignore.case = TRUE)
    if (toupper(chromosome_label) %in% c("X", "Y")) {
        return(list(status = "excluded", reason = "non_autosomal_vmr"))
    }
    chromosome_dir <- paste0("chr_", chromosome_label)
    stem <- paste0(task$start, "_", task$end)
    bed <- file.path(
        vmr_run_dir, "plink_format", chromosome_dir,
        paste0("TOPMed_LIBD-", cohort, ".", stem, ".bed")
    )
    no_snp <- sub("\\.bed$", ".no-snps", bed)
    if (!file.exists(bed)) {
        if (file.exists(no_snp)) {
            return(list(status = "qc_failed",
                        reason = "no_snp_in_prespecified_cis_window"))
        }
        stop("Missing PLINK BED: ", bed)
    }
    prefix <- if (identical(cohort, "AA")) "TOPMed_LIBD.AA" else "TOPMed_LIBD"
    phenotype_path <- file.path(
        vmr_run_dir, "vmr", "phenotypes",
        paste0(task$chrom, "_", stem, "_meth.phen")
    )
    covar_path <- file.path(vmr_run_dir, "covs", chromosome_dir,
                            paste0(prefix, ".covar"))
    qcovar_path <- file.path(vmr_run_dir, "covs", chromosome_dir,
                             paste0(prefix, ".qcovar"))
    for (path in c(phenotype_path, covar_path, qcovar_path)) {
        if (!file.exists(path)) stop("Missing observed input: ", path)
    }

    if (!requireNamespace("bigsnpr", quietly = TRUE)) {
        stop("bigsnpr is required")
    }
    backing <- tempfile(pattern = paste0(backing_tag, "-"))
    rds <- bigsnpr::snp_readBed(bed, backingfile = backing)
    obj <- bigsnpr::snp_attach(rds)
    on.exit(unlink(c(obj$genotypes$backingfile, obj$genotypes$rds),
                   force = TRUE), add = TRUE)
    genotype <- as.matrix(obj$genotypes[])
    map <- obj$map
    if (!all(c("chromosome", "physical.pos") %in% names(map))) {
        stop("PLINK map lacks chromosome or physical.pos")
    }
    window_bp <- 500000L
    window_start <- max(1L, as.integer(task$start) - window_bp)
    window_end <- as.integer(task$end) + window_bp
    map_chr <- sub("^chr", "", as.character(map$chromosome), ignore.case = TRUE)
    in_window <- map_chr == chromosome_label &
        map$physical.pos >= window_start & map$physical.pos <= window_end
    if (!any(in_window)) {
        return(list(status = "qc_failed",
                    reason = "no_snp_in_prespecified_cis_window"))
    }
    genotype <- genotype[, in_window, drop = FALSE]
    snps_in_window <- ncol(genotype)
    missingness <- colMeans(is.na(genotype))
    af <- colMeans(genotype, na.rm = TRUE) / 2
    maf <- pmin(af, 1 - af)
    keep_snp <- is.finite(maf) & maf >= 0.05 &
        is.finite(missingness) & missingness <= 0.05
    ## Gate on the QC'd count regardless, so both callers see the same loci.
    if (sum(keep_snp) < min_cis_variants) {
        return(list(status = "qc_failed",
                    reason = "fewer_than_min_cis_variants",
                    snps_in_window = snps_in_window))
    }
    if (apply_snp_qc) genotype <- genotype[, keep_snp, drop = FALSE]

    fam <- obj$fam[, 1:2, drop = FALSE]
    names(fam) <- c("FID", "IID")
    fam$FID <- as.character(fam$FID)
    fam$IID <- as.character(fam$IID)
    phenotype <- read.table(phenotype_path, header = FALSE,
                            stringsAsFactors = FALSE)
    names(phenotype) <- c("FID", "IID", "phenotype")
    covar <- read.table(covar_path, header = FALSE, stringsAsFactors = FALSE)
    names(covar) <- c("FID", "IID", "sex", "diagnosis")
    qcovar <- read.table(qcovar_path, header = FALSE, stringsAsFactors = FALSE)
    names(qcovar) <- c("FID", "IID", "age")
    metadata <- Reduce(
        function(x, y) merge(x, y, by = c("FID", "IID"), all = FALSE),
        list(fam, phenotype, covar, qcovar)
    )
    metadata$key <- paste(metadata$FID, metadata$IID, sep = "::")
    fam$key <- paste(fam$FID, fam$IID, sep = "::")
    row_index <- match(metadata$key, fam$key)
    if (anyNA(row_index) || anyDuplicated(metadata$key)) {
        stop("Donor alignment failed or produced duplicate IDs")
    }
    genotype <- genotype[row_index, , drop = FALSE]
    keep_sample <- is.finite(as.numeric(metadata$phenotype)) &
        is.finite(as.numeric(metadata$age)) &
        !is.na(metadata$sex) & !is.na(metadata$diagnosis)
    metadata <- metadata[keep_sample, , drop = FALSE]
    genotype <- genotype[keep_sample, , drop = FALSE]
    if (is.finite(expected_n) && nrow(genotype) != as.integer(expected_n)) {
        stop("Observed donor count differs from locked design: ", nrow(genotype),
             " versus ", expected_n)
    }
    y <- as.numeric(metadata$phenotype)
    covariates <- stats::model.matrix(
        ~ age + factor(sex) + factor(diagnosis), data = metadata
    )[, -1L, drop = FALSE]

    list(
        status = "ok", reason = NA_character_,
        genotype = genotype, y = y, covariates = covariates,
        metadata = metadata, snps_in_window = snps_in_window,
        plink_source = normalizePath(bed),
        phenotype_source = normalizePath(phenotype_path)
    )
}

## Replace a locus's observed phenotype with one simulated at a known true PVE,
## holding the real genotype fixed. This is the only difference between the
## observed-regime grid and production Stage 01.
simulate_phenotype_on_observed_genotype <- function(genotype, covariates, h2,
                                                    architecture) {
    n <- nrow(genotype)
    p <- ncol(genotype)
    dosage <- genotype
    if (anyNA(dosage)) {
        means <- colMeans(dosage, na.rm = TRUE)
        means[!is.finite(means)] <- 0
        for (j in seq_len(p)) {
            miss <- is.na(dosage[, j])
            if (any(miss)) dosage[miss, j] <- means[[j]]
        }
    }
    number_causal <- if (architecture == "sparse") {
        min(5L, p)
    } else if (architecture == "oligogenic") {
        min(p, max(10L, ceiling(p * 0.01)))
    } else if (architecture == "polygenic") {
        p
    } else {
        stop("Unknown architecture: ", architecture)
    }
    causal_index <- if (h2 > 0) sample(seq_len(p), number_causal) else integer()
    beta <- numeric(p)
    if (length(causal_index)) {
        beta[causal_index] <- stats::rnorm(length(causal_index))
    }
    genetic_value <- drop(dosage %*% beta)
    if (h2 > 0 && stats::var(genetic_value) > 0) {
        genetic_value <- as.numeric(scale(genetic_value)) * sqrt(h2)
    } else {
        genetic_value[] <- 0
    }
    noise <- stats::rnorm(n)
    if (stats::var(genetic_value) > 0) {
        noise <- stats::residuals(stats::lm(noise ~ genetic_value))
    }
    noise <- as.numeric(scale(noise)) * sqrt(1 - h2)
    residual_phenotype <- genetic_value + noise
    covariate_beta <- seq(0.15, 0.05, length.out = ncol(covariates))
    list(
        phenotype = residual_phenotype + drop(covariates %*% covariate_beta),
        genetic_value = genetic_value,
        causal_index = causal_index,
        realized_h2 = safe_ratio(stats::var(genetic_value),
                                 stats::var(residual_phenotype))
    )
}
