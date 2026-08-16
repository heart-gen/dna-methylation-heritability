#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    region = Sys.getenv("REGION", unset = ""),
    population = Sys.getenv("POPULATION", unset = ""),
    task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = ""),
    repo_root = Sys.getenv(
        "CAL_H2_REPO_ROOT",
        unset = normalizePath(file.path(dirname(script_path), "..", ".."))
    ),
    ## v2: point the adapter at an accepted 01_vmr_catalog run instead of the
    ## legacy vmr-analysis/all_individuals tree. When set, vmr_run_dir supplies
    ## vmr.bed, the per-VMR phenotypes, the covariates, and plink_format, so the
    ## VMR set and the donor set come from one immutable, checksummed source.
    ## Legacy behaviour is preserved when it is empty.
    vmr_run_dir = Sys.getenv("CAL_H2_VMR_RUN_DIR", unset = ""),
    ## Cohort arm: AA (primary) or all_individuals (sensitivity).
    cohort = Sys.getenv("CAL_H2_COHORT", unset = ""),
    plink_root = Sys.getenv("CAL_H2_PLINK_ROOT", unset = ""),
    recovered_plink_root = Sys.getenv("CAL_H2_RECOVERED_PLINK_ROOT", unset = ""),
    phenotype_root = Sys.getenv("CAL_H2_PHENOTYPE_ROOT", unset = ""),
    calibration_model = file.path(dirname(script_path), "..", "_m", "calibration", "elastic-net-calibration.rds"),
    output_root = file.path(dirname(script_path), "..", "_m", "observed"),
    cis_window_bp = "500000",
    maf_min = "0.05",
    snp_missingness_max = "0.05",
    ## config/thresholds.yml cis.min_cis_variants. Kept as a CLI default rather
    ## than read here so this script stays runnable standalone; step_5 passes
    ## the configured value.
    min_cis_variants = "100",
    write_diagnostics = "FALSE",
    include_sex_chromosomes = "FALSE",
    allow_estimator_mismatch = "FALSE"
))
if (!nzchar(cli$region) || !nzchar(cli$population)) {
    stop("region and population are required")
}
task_id <- as_int(cli$task_id, "task_id")
allow_mismatch <- as_bool(cli$allow_estimator_mismatch, "allow_estimator_mismatch")
write_diagnostics <- as_bool(cli$write_diagnostics, "write_diagnostics")

region <- tolower(cli$region)
population <- cli$population

## ---------------------------------------------------------------------------
## v2 input source.
##
## When vmr_run_dir is set, every VMR-derived input comes from that one accepted
## 01_vmr_catalog run: vmr.bed, the per-VMR phenotypes, the covariates, and the
## cis genotype windows. The legacy path below reads from
## vmr-analysis/all_individuals/, whose VMR sets are invalidated by defect V1
## (donor rows misaligned against PC design rows) and are therefore usable only
## for old-versus-new comparison, never as production input.
## ---------------------------------------------------------------------------
using_v2_catalog <- nzchar(cli$vmr_run_dir)
if (using_v2_catalog) {
    if (!dir.exists(cli$vmr_run_dir)) {
        stop("vmr_run_dir does not exist: ", cli$vmr_run_dir)
    }
    v2_manifest <- file.path(cli$vmr_run_dir, "manifest.tsv")
    if (!file.exists(v2_manifest)) {
        stop("vmr_run_dir has no manifest.tsv, so it is not an accepted run: ",
             cli$vmr_run_dir)
    }
    base_dir <- cli$vmr_run_dir
    vmr_file <- file.path(base_dir, "vmr", "vmr.bed")

    ## Carry the upstream identity forward. AGENTS.md 9 requires every output to
    ## record its upstream run IDs and vmr_set_id; AGENTS.md 6 forbids consuming
    ## an upstream result that has not recorded a passing acceptance gate.
    v2_man <- read.table(v2_manifest, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, quote = "", comment.char = "")
    v2_field <- function(f, default = NA_character_) {
        i <- match(f, v2_man$field)
        if (is.na(i)) default else v2_man$value[i]
    }
    upstream_run_id <- v2_field("run_id")
    upstream_vmr_set_id <- v2_field("vmr_set_id")
    upstream_cohort <- v2_field("cohort")
    upstream_region <- v2_field("region")

    if (!nzchar(cli$cohort)) {
        cli$cohort <- upstream_cohort
    } else if (!identical(cli$cohort, upstream_cohort)) {
        stop("cohort mismatch: requested '", cli$cohort, "' but ",
             "vmr_run_dir was built for '", upstream_cohort, "'")
    }
    if (!identical(tolower(upstream_region), region)) {
        stop("region mismatch: requested '", region, "' but vmr_run_dir was ",
             "built for '", upstream_region, "'")
    }

    ## Covariate filename prefix is per-arm (TOPMed_LIBD.AA vs TOPMed_LIBD).
    v2_covar_prefix <- if (identical(cli$cohort, "AA")) {
        "TOPMed_LIBD.AA"
    } else {
        "TOPMed_LIBD"
    }

    message("[v2] upstream run ", upstream_run_id,
            " | vmr_set_id ", upstream_vmr_set_id,
            " | cohort ", cli$cohort, " | region ", region)
} else {
    upstream_run_id <- NA_character_
    upstream_vmr_set_id <- NA_character_
    v2_covar_prefix <- NA_character_
    base_dir <- file.path(cli$repo_root, "vmr-analysis", "all_individuals",
                          region, "_m")
    vmr_file <- file.path(base_dir, "vmr.bed")
}

if (nzchar(cli$plink_root)) {
    plink_base_dir <- file.path(cli$plink_root, region, "_m", "plink_format")
} else {
    local_plink <- file.path(base_dir, "plink_format")
    if (using_v2_catalog) {
        ## No cross-repo fallback: a v2 run must read the genotype windows that
        ## were extracted for its own VMR set, or fail.
        if (!dir.exists(local_plink)) {
            stop("vmr_run_dir has no plink_format/: ", local_plink,
                 "\n  Run 01_vmr_catalog/_h/step_4.sh for this run first.")
        }
        plink_base_dir <- local_plink
    } else {
        alexis_plink <- file.path(
            "/projects/b1213/users/alexis/projects/dna-methylation-heritability",
            "vmr-analysis", "all_individuals", region, "_m", "plink_format"
        )
        plink_base_dir <- if (dir.exists(local_plink)) local_plink else alexis_plink
    }
}
vmrs <- read.table(vmr_file, header = FALSE, stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(vmrs)) stop("task_id is outside vmr.bed")
chromosome <- as.character(vmrs[task_id, 1L])
start <- as.integer(vmrs[task_id, 2L])
end <- as.integer(vmrs[task_id, 3L])
chromosome_label <- sub("^chr", "", chromosome, ignore.case = TRUE)
chromosome_dir <- paste0("chr_", chromosome_label)
stem <- paste0(start, "_", end)
vmr_record <- data.frame(
    task_id = task_id,
    region = region,
    population = population,
    chromosome = chromosome,
    start = start,
    end = end,
    vmr_id = paste(chromosome, start, end, sep = ":"),
    ## Upstream identity travels with every row, so a downstream module can tell
    ## which VMR catalog a number came from without consulting a side file.
    upstream_vmr_run_id = upstream_run_id,
    vmr_set_id = upstream_vmr_set_id,
    cohort = if (nzchar(cli$cohort)) cli$cohort else NA_character_,
    stringsAsFactors = FALSE
)
write_exclusion <- function(reason, source_log = NA_character_) {
    excluded <- cbind(
        vmr_record,
        data.frame(
            exclusion_reason = reason,
            source_log = source_log,
            stringsAsFactors = FALSE
        )
    )
    write_tsv(
        excluded,
        file.path(cli$output_root, region, population, "excluded",
                  sprintf("vmr-%07d.tsv", task_id))
    )
    cat("Excluded", excluded$vmr_id, "because", reason, "\n")
}
write_qc_failure <- function(reason, snps_in_window, snps_after_qc) {
    failure <- cbind(
        vmr_record,
        data.frame(
            qc_failure_reason = reason,
            snps_in_window = snps_in_window,
            snps_after_qc = snps_after_qc,
            stringsAsFactors = FALSE
        )
    )
    write_tsv(
        failure,
        file.path(cli$output_root, region, population, "qc_failures",
                  sprintf("vmr-%07d.tsv", task_id))
    )
    cat("Recorded QC failure for", failure$vmr_id, ":", reason, "\n")
}
if (toupper(chromosome_label) %in% c("X", "Y") &&
    !as_bool(cli$include_sex_chromosomes, "include_sex_chromosomes")) {
    write_exclusion("non_autosomal_vmr")
    quit(save = "no", status = 0L)
}
if (!requireNamespace("bigsnpr", quietly = TRUE)) stop("The bigsnpr package is required")

bed <- file.path(
    plink_base_dir, chromosome_dir,
    paste0("TOPMed_LIBD-", population, ".", stem, ".bed")
)
if (nzchar(cli$recovered_plink_root)) {
    recovered_bed <- file.path(
        cli$recovered_plink_root, region, "_m", "plink_format",
        chromosome_dir,
        paste0("TOPMed_LIBD-", population, ".", stem, ".bed")
    )
    recovered_no_snp <- sub("\\.bed$", ".no-snps", recovered_bed)
    if (file.exists(recovered_bed) || file.exists(recovered_no_snp)) {
        bed <- recovered_bed
    }
}
bed_log <- sub("\\.bed$", ".log", bed)
no_snp_marker <- sub("\\.bed$", ".no-snps", bed)
if (!file.exists(bed)) {
    upstream_bed <- file.path(
        plink_base_dir, chromosome_dir,
        paste0("TOPMed_LIBD-", population, ".", stem, ".bed")
    )
    upstream_log <- sub("\\.bed$", ".log", upstream_bed)
    no_variants <- file.exists(no_snp_marker) ||
        (file.exists(upstream_log) && any(grepl(
            "No variants remaining after main filters",
            readLines(upstream_log, warn = FALSE), fixed = TRUE
        )))
    if (no_variants) {
        write_qc_failure("no_snp_in_prespecified_cis_window", 0L, 0L)
        quit(save = "no", status = 0L)
    }
    stop("Required input is missing: ", bed)
}
if (using_v2_catalog) {
    ## 01_vmr_catalog writes one flat phenotypes/ directory keyed by
    ## chr_start_end, and there is no fallback tree.
    phenotype_file <- file.path(
        base_dir, "vmr", "phenotypes",
        paste0(chromosome, "_", stem, "_meth.phen")
    )
    covar_file <- file.path(base_dir, "covs", chromosome_dir,
                            paste0(v2_covar_prefix, ".covar"))
    qcovar_file <- file.path(base_dir, "covs", chromosome_dir,
                             paste0(v2_covar_prefix, ".qcovar"))
} else {
local_phenotype <- file.path(
    base_dir, "vmr", chromosome_dir, paste0(stem, "_meth.phen")
)
if (nzchar(cli$phenotype_root)) {
    fallback_phenotype <- file.path(
        cli$phenotype_root, region, "_m", "vmr", chromosome_dir,
        paste0(stem, "_meth.phen")
    )
} else {
    fallback_phenotype <- file.path(
        "/projects/b1213/users/alexis/projects/dna-methylation-heritability",
        "vmr-analysis", "all_individuals", region, "_m", "vmr",
        chromosome_dir, paste0(stem, "_meth.phen")
    )
}
phenotype_file <- if (file.exists(local_phenotype)) {
    local_phenotype
} else {
    fallback_phenotype
}
covar_file <- file.path(base_dir, "covs", "TOPMed_LIBD.covar")
qcovar_file <- file.path(base_dir, "covs", "TOPMed_LIBD.qcovar")
}
for (path in c(phenotype_file, covar_file, qcovar_file, cli$calibration_model)) {
    if (!file.exists(path)) stop("Required input is missing: ", path)
}

backing <- tempfile(pattern = paste0("vmr-", task_id, "-"))
rds <- bigsnpr::snp_readBed(bed, backingfile = backing)
big_snp <- bigsnpr::snp_attach(rds)
on.exit({
    unlink(c(big_snp$genotypes$backingfile, big_snp$genotypes$rds), force = TRUE)
}, add = TRUE)
genotype <- as.matrix(big_snp$genotypes[])
map <- big_snp$map
if (!all(c("chromosome", "physical.pos") %in% names(map))) {
    stop("PLINK map does not contain chromosome and physical.pos columns")
}
cis_window_bp <- as_int(cli$cis_window_bp, "cis_window_bp")
window_start <- max(1L, start - cis_window_bp)
window_end <- end + cis_window_bp
map_chromosome <- sub("^chr", "", as.character(map$chromosome), ignore.case = TRUE)
in_window <- map_chromosome == chromosome_label &
    map$physical.pos >= window_start & map$physical.pos <= window_end
if (!any(in_window)) {
    write_qc_failure("no_snp_in_prespecified_cis_window", 0L, 0L)
    quit(save = "no", status = 0L)
}
genotype <- genotype[, in_window, drop = FALSE]
input_snps <- ncol(genotype)
missingness <- colMeans(is.na(genotype))
allele_frequency <- colMeans(genotype, na.rm = TRUE) / 2
maf <- pmin(allele_frequency, 1 - allele_frequency)
maf_min <- as_num(cli$maf_min, "maf_min")
missingness_max <- as_num(cli$snp_missingness_max, "snp_missingness_max")
pass_qc <- is.finite(maf) & maf >= maf_min &
    is.finite(missingness) & missingness <= missingness_max
genotype <- genotype[, pass_qc, drop = FALSE]
if (ncol(genotype) < 2L) {
    write_qc_failure("fewer_than_two_snps_after_qc", input_snps, ncol(genotype))
    quit(save = "no", status = 0L)
}

## Minimum cis variants (config/thresholds.yml cis.min_cis_variants).
##
## A locus with a handful of cis SNPs still produces a number, which is exactly
## the problem: the elastic net fits, the calibration maps it, and nothing
## downstream marks it as resting on almost no predictors. The AA dlpfc catalog
## has 19 such VMRs, six of them with a single variant. That is not the
## estimator the calibration was built against, so the estimate is
## uninterpretable rather than merely imprecise.
##
## Applied AFTER MAF/missingness QC, matching the methods text ("fewer than 100
## after filtering"). Recorded as a QC failure, never a silent drop -- these
## loci cluster by genomic feature, so the analyzed set is not a uniform sample
## of the catalog and the exclusion count has to be reportable.
min_cis_variants <- as_int(cli$min_cis_variants, "min_cis_variants")
if (ncol(genotype) < min_cis_variants) {
    message("[qc] ", ncol(genotype), " cis SNPs after QC (< ", min_cis_variants,
            "); excluding ", chromosome, ":", start, "-", end)
    write_qc_failure("fewer_than_min_cis_variants", input_snps, ncol(genotype))
    quit(save = "no", status = 0L)
}
fam <- big_snp$fam[, 1:2, drop = FALSE]
names(fam) <- c("FID", "IID")
fam$FID <- as.character(fam$FID)
fam$IID <- as.character(fam$IID)

phenotype <- read.table(phenotype_file, header = FALSE, stringsAsFactors = FALSE)
names(phenotype) <- c("FID", "IID", "phenotype")
covar <- read.table(covar_file, header = FALSE, stringsAsFactors = FALSE)
names(covar) <- c("FID", "IID", "sex", "diagnosis")
qcovar <- read.table(qcovar_file, header = FALSE, stringsAsFactors = FALSE)
names(qcovar) <- c("FID", "IID", "age")
metadata <- Reduce(function(x, y) merge(x, y, by = c("FID", "IID"), all = FALSE),
                   list(fam, phenotype, covar, qcovar))
metadata$key <- paste(metadata$FID, metadata$IID, sep = "::")
fam$key <- paste(fam$FID, fam$IID, sep = "::")
row_index <- match(metadata$key, fam$key)
if (anyNA(row_index)) stop("Failed to align phenotype/covariate samples to PLINK FAM")
genotype <- genotype[row_index, , drop = FALSE]
keep_sample <- is.finite(as.numeric(metadata$phenotype)) &
    is.finite(as.numeric(metadata$age)) &
    !is.na(metadata$sex) & !is.na(metadata$diagnosis)
metadata <- metadata[keep_sample, , drop = FALSE]
genotype <- genotype[keep_sample, , drop = FALSE]
covariates <- stats::model.matrix(
    ~ age + factor(sex) + factor(diagnosis), data = metadata
)[, -1L, drop = FALSE]

## AGENTS.md 7.2: "preserve the accepted calibration model and checksum".
## The frozen model is the one that passed all seven acceptance criteria in
## ajhg-calibration-v4-independent-validation. Verifying its SHA-256 here means
## a silently swapped or corrupted model cannot produce h2_en_calibrated values
## that look ordinary but come from an unaccepted calibration.
CALIBRATION_MODEL_SHA256 <-
    "bbe9f9f3e897b19c536078c20e6bd50a2f5ea385ab1c1258039974ced855e389"
observed_model_sha <- {
    out <- suppressWarnings(tryCatch(
        system2("sha256sum", shQuote(normalizePath(cli$calibration_model)),
                stdout = TRUE, stderr = FALSE),
        error = function(e) NA_character_))
    if (length(out) == 0 || is.na(out[[1L]])) NA_character_ else sub(" .*$", "", out[[1L]])
}
if (is.na(observed_model_sha)) {
    stop("Could not checksum the calibration model at ", cli$calibration_model)
}
if (!identical(observed_model_sha, CALIBRATION_MODEL_SHA256)) {
    stop("Calibration model checksum mismatch.\n  expected ",
         CALIBRATION_MODEL_SHA256, "\n  observed ", observed_model_sha,
         "\n  This is not the accepted model from ",
         "ajhg-calibration-v4-independent-validation. Changing the estimator, ",
         "alpha grid, folds, repeats, lambda rule, screen, or raw metric ",
         "requires full recalibration (AGENTS.md 7.2).")
}

model <- readRDS(cli$calibration_model)
settings <- model$estimator_settings
outer_folds <- as.integer(settings$outer_folds)
outer_repeats <- as.integer(settings$outer_repeats)
inner_folds <- as.integer(settings$inner_folds)
alpha_grid <- split_numeric(settings$alpha_grid)
lambda_rule <- as.character(settings$lambda_rule)
max_features <- as.integer(settings$max_features)

requested_settings <- c("outer_folds", "outer_repeats", "inner_folds",
                        "alpha_grid", "lambda_rule", "max_features")
overrides <- intersect(requested_settings, names(cli))
if (length(overrides) && !allow_mismatch) {
    stop("Estimator-setting overrides require --allow-estimator-mismatch=TRUE; recalibration is preferred")
}
if (length(overrides)) {
    if ("outer_folds" %in% overrides) outer_folds <- as_int(cli$outer_folds, "outer_folds")
    if ("outer_repeats" %in% overrides) outer_repeats <- as_int(cli$outer_repeats, "outer_repeats")
    if ("inner_folds" %in% overrides) inner_folds <- as_int(cli$inner_folds, "inner_folds")
    if ("alpha_grid" %in% overrides) alpha_grid <- split_numeric(cli$alpha_grid)
    if ("lambda_rule" %in% overrides) lambda_rule <- cli$lambda_rule
    if ("max_features" %in% overrides) max_features <- as_int(cli$max_features, "max_features")
}

fit <- crossfit_elastic_net(
    genotype = genotype,
    phenotype = as.numeric(metadata$phenotype),
    covariates = covariates,
    outer_folds = outer_folds,
    outer_repeats = outer_repeats,
    inner_folds = inner_folds,
    alpha_grid = alpha_grid,
    lambda_rule = lambda_rule,
    max_features = max_features,
    seed = 20250805L + task_id * 1009L,
    keep_predictions = write_diagnostics
)
fit$metrics$ld_metric <- adjacent_ld_metric(genotype)
he_metrics <- haseman_elston(
    genotype, as.numeric(metadata$phenotype), covariates
)
calibrated <- predict_calibration(model, cbind(fit$metrics, he_metrics))
summary <- cbind(data.frame(
    task_id = task_id,
    region = region,
    population = population,
    chromosome = chromosome,
    start = start,
    end = end,
    vmr_id = paste(chromosome, start, end, sep = ":"),
    samples = nrow(genotype),
    cis_window_bp = cis_window_bp,
    snps_in_window = input_snps,
    snps_after_qc = ncol(genotype),
    maf_min = maf_min,
    snp_missingness_max = missingness_max,
    min_cis_variants = min_cis_variants,
        plink_source = normalizePath(bed),
    phenotype_source = normalizePath(phenotype_file),
    calibration_model = normalizePath(cli$calibration_model),
    stringsAsFactors = FALSE
), fit$metrics, he_metrics, calibrated)

output_dir <- file.path(cli$output_root, region, population)
write_tsv(summary, file.path(output_dir, "summary", sprintf("vmr-%07d.tsv", task_id)))
if (write_diagnostics) {
    write_tsv(fit$folds, file.path(output_dir, "folds", sprintf("vmr-%07d.tsv", task_id)))
    predictions <- cbind(metadata[, c("FID", "IID")], fit$predictions)
    write_tsv(predictions, file.path(output_dir, "predictions", sprintf("vmr-%07d.tsv", task_id)))
}
cat("Completed", summary$vmr_id, "with calibration status", summary$calibration_status, "\n")
