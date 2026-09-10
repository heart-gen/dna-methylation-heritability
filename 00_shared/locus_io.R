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
##
## `cohort` is the CELL token. Two things are derived from it, and they are not
## the same thing once discovery and estimation are separated (AGENTS.md 7.7,
## PI 2026-09-06):
##
##   * the per-VMR BED is named by the ESTIMATION GROUP -- the donors the SNP
##     model is fit in. For a bare arm the group equals the cohort, so every
##     path resolved for AA / all_individuals is unchanged.
##   * the shared covariate files are named by the CATALOG COHORT, because they
##     were written by that arm's Module 01 run. `covar_prefix` overrides this;
##     when NULL it falls back to exactly the pre-2026-09-10 ternary.
##
## MAF and missingness QC below run AFTER the group-restricted BED is loaded, so
## the >= maf_min filter is computed within the estimation group. That is the
## whole point of extracting per-group BEDs rather than subsetting donors here.
##
## If `covs/genotype_pcs.tsv` exists in `vmr_run_dir`, its snpPC columns are
## appended to the covariate matrix. 01b_estimation_cells writes that file from
## a within-group PCA. When it is absent the covariate matrix is byte-identical
## to the pre-2026-09-10 behaviour, so the sealed accepted runs reproduce.
## Read one locus's donor-level phenotype and covariate tables.
##
## Extracted from load_observed_locus() so that callers which model methylation
## WITHOUT genotype -- 10_environmental_exploratory fits
## `meth ~ exposure + covariates` and never touches a dosage -- reuse this
## implementation rather than writing a second one (AGENTS.md 5.3: "Refactor
## shared logic into parameterized functions"). The covariate model is defined
## in exactly one place, so an exposure association and a local-genetic-control
## estimate condition on the same thing by construction.
##
## Returns the tables unmerged, because the two callers join them against
## different donor universes: load_observed_locus() inner-joins against the
## locus BED's `fam`, while a genotype-free caller uses the phenotype file's own
## donors. Merging here would force one of those choices on both.
##
## `covs/genotype_pcs.tsv` is optional. 01b_estimation_cells writes it from a
## within-group PCA; when it is absent the returned covariate set is exactly the
## pre-2026-09-10 age/sex/diagnosis, which is why sealed runs reproduce.
load_locus_phenotype <- function(task, vmr_run_dir, cohort = NULL,
                                 covar_prefix = NULL) {
    chromosome_label <- sub("^chr", "", task$chrom, ignore.case = TRUE)
    chromosome_dir <- paste0("chr_", chromosome_label)
    stem <- paste0(task$start, "_", task$end)
    prefix <- if (!is.null(covar_prefix)) {
        covar_prefix
    } else if (identical(cohort, "AA")) {
        "TOPMed_LIBD.AA"
    } else {
        "TOPMed_LIBD"
    }
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

    phenotype <- read.table(phenotype_path, header = FALSE,
                            stringsAsFactors = FALSE)
    names(phenotype) <- c("FID", "IID", "phenotype")
    covar <- read.table(covar_path, header = FALSE, stringsAsFactors = FALSE)
    names(covar) <- c("FID", "IID", "sex", "diagnosis")
    qcovar <- read.table(qcovar_path, header = FALSE, stringsAsFactors = FALSE)
    names(qcovar) <- c("FID", "IID", "age")

    pc_path <- file.path(vmr_run_dir, "covs", "genotype_pcs.tsv")
    pc_names <- character(0)
    pcs <- NULL
    if (file.exists(pc_path)) {
        pcs <- utils::read.delim(pc_path, header = TRUE, colClasses = "character",
                                 stringsAsFactors = FALSE)
        if (!all(c("FID", "IID") %in% names(pcs))) {
            stop("genotype_pcs.tsv lacks FID/IID: ", pc_path)
        }
        pc_names <- grep("^snpPC[0-9]+$", names(pcs), value = TRUE)
        if (!length(pc_names)) {
            stop("genotype_pcs.tsv carries no snpPC columns: ", pc_path)
        }
        ## Order by index, not by the file's column order, so the covariate
        ## matrix is reproducible whatever wrote the file.
        pc_names <- pc_names[order(as.integer(sub("^snpPC", "", pc_names)))]
        pcs <- pcs[, c("FID", "IID", pc_names), drop = FALSE]
        for (nm in pc_names) pcs[[nm]] <- as.numeric(pcs[[nm]])
        if (anyNA(pcs[, pc_names, drop = FALSE])) {
            stop("Non-numeric or missing genotype PC values in ", pc_path)
        }
        if (anyDuplicated(paste(pcs$FID, pcs$IID, sep = "::"))) {
            stop("Duplicate donors in ", pc_path)
        }
    }

    list(phenotype = phenotype, covar = covar, qcovar = qcovar,
         pcs = pcs, pc_names = pc_names,
         pc_path = if (length(pc_names)) pc_path else NA_character_,
         phenotype_source = normalizePath(phenotype_path),
         covar_source = normalizePath(covar_path))
}

load_observed_locus <- function(task, cohort, vmr_run_dir, min_cis_variants,
                                expected_n = NA_integer_,
                                backing_tag = "lgv",
                                apply_snp_qc = TRUE,
                                covar_prefix = NULL,
                                estimation_group = NULL) {
    ## Resolve the cell without requiring callers to have done it. Config is
    ## unavailable in some unit-test harnesses, so fall back to the identity
    ## mapping (group == cohort), which is the pre-2026-09-10 behaviour.
    if (is.null(estimation_group)) {
        estimation_group <- tryCatch(
            parse_cell(cohort)$estimation_group,
            error = function(e) cohort
        )
    }
    chromosome_label <- sub("^chr", "", task$chrom, ignore.case = TRUE)
    if (toupper(chromosome_label) %in% c("X", "Y")) {
        return(list(status = "excluded", reason = "non_autosomal_vmr"))
    }
    chromosome_dir <- paste0("chr_", chromosome_label)
    stem <- paste0(task$start, "_", task$end)
    bed <- file.path(
        vmr_run_dir, "plink_format", chromosome_dir,
        paste0("TOPMed_LIBD-", estimation_group, ".", stem, ".bed")
    )
    no_snp <- sub("\\.bed$", ".no-snps", bed)
    if (!file.exists(bed)) {
        if (file.exists(no_snp)) {
            return(list(status = "qc_failed",
                        reason = "no_snp_in_prespecified_cis_window"))
        }
        stop("Missing PLINK BED: ", bed)
    }
    pheno <- load_locus_phenotype(task = task, vmr_run_dir = vmr_run_dir,
                                  cohort = cohort, covar_prefix = covar_prefix)

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
    ## A donor genotyped at this locus but absent from the PC table would be
    ## dropped by the inner merge below without a word, changing n silently.
    ## That is the V1 class of defect; make it a stop. It needs `fam`, so it
    ## stays here rather than moving into load_locus_phenotype(), which has no
    ## genotype to compare against.
    tables <- list(fam, pheno$phenotype, pheno$covar, pheno$qcovar)
    pc_names <- pheno$pc_names
    if (!is.null(pheno$pcs)) {
        fam_key <- paste(fam$FID, fam$IID, sep = "::")
        pc_key <- paste(pheno$pcs$FID, pheno$pcs$IID, sep = "::")
        missing_pc <- setdiff(fam_key, pc_key)
        if (length(missing_pc)) {
            stop("Donors present in the locus BED but absent from ",
                 pheno$pc_path, ": ",
                 paste(head(missing_pc, 10), collapse = ", "),
                 if (length(missing_pc) > 10)
                     paste0(" (and ", length(missing_pc) - 10, " more)"))
        }
        tables <- c(tables, list(pheno$pcs))
    }

    metadata <- Reduce(
        function(x, y) merge(x, y, by = c("FID", "IID"), all = FALSE),
        tables
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
    ## Build the formula rather than branching on two literals: with no PC file
    ## this is exactly `~ age + factor(sex) + factor(diagnosis)`, the
    ## pre-2026-09-10 model, term for term and in the same order.
    terms <- c("age", "factor(sex)", "factor(diagnosis)", pc_names)
    covariates <- stats::model.matrix(
        stats::reformulate(terms), data = metadata
    )[, -1L, drop = FALSE]

    list(
        status = "ok", reason = NA_character_,
        genotype = genotype, y = y, covariates = covariates,
        metadata = metadata, snps_in_window = snps_in_window,
        plink_source = normalizePath(bed),
        phenotype_source = pheno$phenotype_source
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
