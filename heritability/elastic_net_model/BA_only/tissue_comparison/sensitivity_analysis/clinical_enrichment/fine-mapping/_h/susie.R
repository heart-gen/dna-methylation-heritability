#!/usr/bin/env Rscript

# =====================================================
# SuSiE Fine-Mapping Script (Summary Statistics + LD)
# - Automatically detects case-control N
# - Uses effective sample size when appropriate
# - Saves PIPs and credible sets
# =====================================================

suppressPackageStartupMessages({
  library(susieR)
  library(data.table)
})

# ----------------------------
# Parse arguments
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript run_susie_finemap.R <sumstats.txt> <ld_matrix.txt> <output_prefix> [default_N]")
}

SUMSTATS_FILE <- args[1]
LD_FILE       <- args[2]
OUT_PREFIX    <- args[3]
DEFAULT_N     <- ifelse(length(args) >= 4, as.numeric(args[4]), NA)

cat("====================================\n")
cat("Running SuSiE fine-mapping\n")
cat("Sumstats:", SUMSTATS_FILE, "\n")
cat("LD matrix:", LD_FILE, "\n")
cat("====================================\n")

# ----------------------------
# Load summary statistics
# ----------------------------
sumstats <- fread(SUMSTATS_FILE)

required_cols <- c("rsid", "beta", "standard_error")
missing_cols <- setdiff(required_cols, colnames(sumstats))

if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Rename to standardized names
setnames(sumstats,
         old = c("rsid", "beta", "standard_error"),
         new = c("SNP", "BETA", "SE"))

# Compute Z-scores
sumstats[, Z := BETA / SE]

# Remove SNPs with missing or infinite Z
sumstats <- sumstats[is.finite(Z)]

# ----------------------------
# Determine Sample Size
# ----------------------------

if (all(c("N_cases", "N_controls") %in% colnames(sumstats))) {

  # Compute per-SNP effective sample size
  sumstats[, N_eff := 4 / ((1 / N_cases) + (1 / N_controls))]

  # Remove invalid values
  sumstats <- sumstats[is.finite(N_eff) & N_eff > 1]

  # Use median effective N across SNPs
  N_scalar <- median(sumstats$N_eff, na.rm = TRUE)

  cat("Using median effective sample size:", N_scalar, "\n")

} else if ("N" %in% colnames(sumstats)) {

  N_scalar <- median(sumstats$N, na.rm = TRUE)
  cat("Using median N:", N_scalar, "\n")

} else if (!is.na(DEFAULT_N)) {

  N_scalar <- DEFAULT_N
  cat("Using provided default N:", N_scalar, "\n")

} else {

  stop("No valid sample size information found.")

}

if (N_scalar <= 1) {
  stop("Sample size must be > 1")
}

# ----------------------------
# Load LD matrix
# ----------------------------
R <- as.matrix(fread(LD_FILE, header = FALSE))

if (nrow(R) != nrow(sumstats)) {
  stop("LD matrix and summary stats have different number of SNPs!")
}

# Identify SNPs with NA in LD
na_snps <- apply(R, 1, function(x) any(is.na(x)))

if (any(na_snps)) {
  cat("Removing", sum(na_snps), "SNPs with NA values in LD matrix\n")
  R <- R[!na_snps, !na_snps]
  sumstats <- sumstats[!na_snps]
}

# Force symmetry
R <- (R + t(R)) / 2

# Regularize slightly
diag(R) <- diag(R) + 1e-8

# ----------------------------
# Run SuSiE
# ----------------------------
cat("Running susie_rss...\n")

susie_fit <- susie_rss(
  z = sumstats$Z,
  R = R,
  n = N_scalar,
  L = 10,
  estimate_residual_variance = TRUE
)

cat("SuSiE finished.\n")

# ----------------------------
# Save Posterior Inclusion Probabilities
# ----------------------------
results <- data.table(
  SNP = sumstats$SNP,
  PIP = susie_fit$pip
)

fwrite(results, paste0(OUT_PREFIX, ".pip.txt"), sep = "\t")

# ----------------------------
# Save Credible Sets
# ----------------------------
cs <- susie_get_cs(susie_fit)

if (!is.null(cs$cs)) {
  cs_out <- data.table(
    CS = rep(seq_along(cs$cs), lengths(cs$cs)),
    SNP = sumstats$SNP[unlist(cs$cs)]
  )
  fwrite(cs_out, paste0(OUT_PREFIX, ".credible_sets.txt"), sep = "\t")
} else {
  cat("No credible sets detected.\n")
}

# ----------------------------
# Save Full Model
# ----------------------------
saveRDS(susie_fit, paste0(OUT_PREFIX, ".susie_model.rds"))

cat("====================================\n")
cat("Outputs written:\n")
cat(paste0(OUT_PREFIX, ".pip.txt\n"))
cat(paste0(OUT_PREFIX, ".credible_sets.txt\n"))
cat(paste0(OUT_PREFIX, ".susie_model.rds\n"))
cat("Done.\n")
cat("====================================\n")