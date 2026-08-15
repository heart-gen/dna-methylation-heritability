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
    plink_root = Sys.getenv("CAL_H2_PLINK_ROOT", unset = ""),
    recovered_plink_root = Sys.getenv("CAL_H2_RECOVERED_PLINK_ROOT", unset = ""),
    phenotype_root = Sys.getenv("CAL_H2_PHENOTYPE_ROOT", unset = ""),
    calibration_model = file.path(dirname(script_path), "..", "_m", "calibration", "elastic-net-calibration.rds"),
    output_root = file.path(dirname(script_path), "..", "_m", "observed"),
    cis_window_bp = "500000",
    maf_min = "0.05",
    snp_missingness_max = "0.05",
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
base_dir <- file.path(cli$repo_root, "vmr-analysis", "all_individuals", region, "_m")
if (nzchar(cli$plink_root)) {
    plink_base_dir <- file.path(cli$plink_root, region, "_m", "plink_format")
} else {
    local_plink <- file.path(base_dir, "plink_format")
    alexis_plink <- file.path(
        "/projects/b1213/users/alexis/projects/dna-methylation-heritability",
        "vmr-analysis", "all_individuals", region, "_m", "plink_format"
    )
    plink_base_dir <- if (dir.exists(local_plink)) local_plink else alexis_plink
}
vmr_file <- file.path(base_dir, "vmr.bed")
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
