#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    input = file.path(dirname(script_path), "..", "_m", "combined",
                      "calibrated-local-h2-all-cells.tsv"),
    output = file.path(dirname(script_path), "..", "_m", "audit",
                       "high-signal-audit-manifest.tsv"),
    n_loci_per_cell = "10",
    seed_repeats = "5",
    base_seed = "20260817"
))

observed <- read_tsv(cli$input)
required <- c("task_id", "region", "population", "rho2_oof", "r2_oof",
              "positive_signal", "calibration_status", "upstream_vmr_run_id",
              "vmr_set_id")
missing <- setdiff(required, names(observed))
if (length(missing)) stop("Observed table lacks: ", paste(missing, collapse = ", "))

n_loci <- as_int(cli$n_loci_per_cell, "n_loci_per_cell")
seed_repeats <- as_int(cli$seed_repeats, "seed_repeats")
base_seed <- as_int(cli$base_seed, "base_seed")
high <- observed[
    observed$calibration_status == "raw_metric_extrapolation" &
        observed$positive_signal %in% TRUE &
        is.finite(observed$rho2_oof) & is.finite(observed$r2_oof),
    , drop = FALSE
]
if (!nrow(high)) stop("No high-signal raw-metric extrapolations were found")

cell_key <- interaction(high$population, high$region, drop = TRUE,
                        lex.order = TRUE)
selected <- do.call(rbind, lapply(split(high, cell_key), function(x) {
    x <- x[order(x$rho2_oof, x$task_id), , drop = FALSE]
    take <- unique(round(seq(1, nrow(x), length.out = min(n_loci, nrow(x)))))
    x[take, , drop = FALSE]
}))
rownames(selected) <- NULL

records <- list()
record_id <- 1L
for (i in seq_len(nrow(selected))) {
    locus <- selected[i, , drop = FALSE]
    sensitivities <- "baseline"
    if (locus$region %in% c("dlpfc", "hippocampus")) {
        sensitivities <- c(sensitivities, "exclude_Br1105")
    }
    population_offset <- if (identical(locus$population, "AA")) 0L else 1000003L
    region_offset <- match(locus$region, c("caudate", "dlpfc", "hippocampus")) * 100003L
    for (sensitivity in sensitivities) {
        for (seed_repeat in seq_len(seed_repeats)) {
            records[[record_id]] <- data.frame(
                audit_id = record_id,
                population = locus$population,
                region = locus$region,
                task_id = locus$task_id,
                vmr_id = locus$vmr_id,
                upstream_vmr_run_id = locus$upstream_vmr_run_id,
                vmr_set_id = locus$vmr_set_id,
                original_r2_oof = locus$r2_oof,
                original_rho2_oof = locus$rho2_oof,
                sensitivity = sensitivity,
                excluded_fids = if (sensitivity == "exclude_Br1105") {
                    "Br1105"
                } else {
                    # A non-empty sentinel is required because Bash treats tab
                    # as IFS whitespace and collapses an empty field, shifting
                    # seed_repeat and fold_seed one column to the left.
                    "NONE"
                },
                seed_repeat = seed_repeat,
                fold_seed = base_seed + population_offset + region_offset +
                    as.integer(locus$task_id) * 1009L + seed_repeat * 9173L,
                stringsAsFactors = FALSE
            )
            record_id <- record_id + 1L
        }
    }
}
manifest <- do.call(rbind, records)
if (anyDuplicated(manifest$audit_id)) stop("Duplicate audit IDs")
write_tsv(manifest, cli$output)
cat("Selected", nrow(selected), "high-signal loci and wrote", nrow(manifest),
    "fold/relatedness audit tasks to", normalizePath(cli$output), "\n")
