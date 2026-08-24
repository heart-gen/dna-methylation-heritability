#!/usr/bin/env Rscript
## Evaluation-only manifest on the locked Module 02 design grid.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))
source(file.path(dirname(script_path), "bslmm_pilot_functions.R"))

cli <- parse_cli(list(
    config = file.path(dirname(script_path), "..", "config", "bslmm-validation.tsv"),
    output = file.path(dirname(script_path), "..", "_m", "config",
                       "bslmm-validation-scenarios.tsv")
))
settings <- read_pilot_settings(cli$config)
required <- c(
    "sample_sizes", "num_snps", "ld_rho", "architectures", "h2_values",
    "evaluation_replicates_per_stratum", "base_seed", "seed_offset",
    "bslmm_burn_in", "bslmm_sampling", "bslmm_rpace", "bslmm_mode", "gemma_bin"
)
missing <- setdiff(required, names(settings))
if (length(missing)) stop("Missing settings: ", paste(missing, collapse = ", "))

architectures <- trimws(strsplit(settings$architectures, ",", fixed = TRUE)[[1L]])
h2_values <- split_numeric(settings$h2_values)
design_cells <- expand.grid(
    architecture = factor(architectures, levels = architectures),
    true_h2 = h2_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
design_cells$architecture <- as.character(design_cells$architecture)
reps <- as_int(settings$evaluation_replicates_per_stratum,
               "evaluation_replicates_per_stratum")
if (reps < nrow(design_cells) || reps %% nrow(design_cells) != 0L) {
    stop(
        "evaluation_replicates_per_stratum must be a positive multiple of ",
        nrow(design_cells), " (architecture-by-h2 cells)"
    )
}

strata <- expand.grid(
    n = as.integer(split_numeric(settings$sample_sizes)),
    num_snps = as.integer(split_numeric(settings$num_snps)),
    ld_rho = split_numeric(settings$ld_rho),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
base_seed <- as_int(settings$base_seed, "base_seed") +
    as_int(settings$seed_offset, "seed_offset")

records <- vector("list", nrow(strata))
for (i in seq_len(nrow(strata))) {
    set.seed(base_seed + i)
    cell_index <- rep(seq_len(nrow(design_cells)), length.out = reps)
    cell_index <- sample(cell_index, length(cell_index), replace = FALSE)
    cell_occurrence <- ave(cell_index, cell_index, FUN = seq_along)
    records[[i]] <- data.frame(
        split = "evaluation",
        n = strata$n[[i]],
        num_snps = strata$num_snps[[i]],
        ld_rho = strata$ld_rho[[i]],
        architecture = design_cells$architecture[cell_index],
        true_h2 = design_cells$true_h2[cell_index],
        replicate = cell_occurrence,
        stringsAsFactors = FALSE
    )
}
manifest <- do.call(rbind, records)
manifest$scenario_id <- seq_len(nrow(manifest))
manifest$seed <- base_seed + manifest$scenario_id * 1009L
manifest$bslmm_mode <- as_int(settings$bslmm_mode, "bslmm_mode")
manifest$bslmm_burn_in <- as_int(settings$bslmm_burn_in, "bslmm_burn_in")
manifest$bslmm_sampling <- as_int(settings$bslmm_sampling, "bslmm_sampling")
manifest$bslmm_rpace <- as_int(settings$bslmm_rpace, "bslmm_rpace")
manifest$gemma_bin <- settings$gemma_bin
manifest$positive_signal_rule <- settings$positive_signal_rule %||% "mcmc_lower_gt_0"
manifest <- manifest[, c(
    "scenario_id", "split", "seed", "n", "num_snps", "ld_rho", "architecture",
    "true_h2", "replicate", "bslmm_mode", "bslmm_burn_in", "bslmm_sampling",
    "bslmm_rpace", "gemma_bin", "positive_signal_rule"
)]
write_tsv(manifest, cli$output)
cat(
    "Wrote", nrow(manifest), "BSLMM validation scenarios (",
    length(unique(paste(manifest$n, manifest$num_snps, manifest$ld_rho))),
    "N/SNP/LD strata) to", normalizePath(cli$output), "\n"
)
