#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "00_functions.R"))

cli <- parse_cli(list(
    manifest = file.path(dirname(script_path), "..", "_m", "audit",
                         "high-signal-audit-manifest.tsv"),
    summaries = file.path(dirname(script_path), "..", "_m", "audit",
                          "summaries"),
    output_dir = file.path(dirname(script_path), "..", "_m", "audit", "summary")
))

manifest <- read_tsv(cli$manifest)
files <- list.files(cli$summaries, pattern = "^audit-[0-9]+\\.tsv$",
                    full.names = TRUE)
if (length(files) != nrow(manifest)) {
    stop("Expected ", nrow(manifest), " summaries but found ", length(files))
}
observed <- do.call(rbind, lapply(files, read_tsv))
audit_ids <- as.integer(sub("^audit-0*([0-9]+)\\.tsv$", "\\1", basename(files)))
observed$audit_id <- audit_ids
audit <- merge(manifest, observed, by = "audit_id", suffixes = c("_planned", ""),
               all.x = TRUE, sort = FALSE)
if (anyNA(audit$rho2_oof)) stop("At least one audit task lacks a finite rho2_oof")
if (any(as.integer(audit$fold_seed_planned) != as.integer(audit$fold_seed))) {
    stop("At least one audit task did not use its planned fold seed")
}
expected_exclusions <- ifelse(
    audit$excluded_fids_planned == "NONE", NA_character_,
    audit$excluded_fids_planned
)
actual_exclusions <- audit$excluded_fids
actual_exclusions[is.na(actual_exclusions) | !nzchar(actual_exclusions) |
                  actual_exclusions == "NA"] <- NA_character_
if (any(xor(is.na(expected_exclusions), is.na(actual_exclusions))) ||
        any(expected_exclusions[!is.na(expected_exclusions)] !=
            actual_exclusions[!is.na(expected_exclusions)])) {
    stop("At least one audit task did not use its planned donor exclusion")
}

baseline <- audit[audit$sensitivity == "baseline", , drop = FALSE]
locus_key <- interaction(baseline$population_planned, baseline$region_planned,
                         baseline$task_id_planned, drop = TRUE, lex.order = TRUE)
locus_summary <- do.call(rbind, lapply(split(baseline, locus_key), function(x) {
    data.frame(
        population = x$population_planned[[1L]],
        region = x$region_planned[[1L]],
        task_id = x$task_id_planned[[1L]],
        vmr_id = x$vmr_id_planned[[1L]],
        seed_repeats = nrow(x),
        positive_r2_fraction = mean(x$r2_oof > 0),
        rho2_median = median(x$rho2_oof),
        rho2_min = min(x$rho2_oof),
        rho2_max = max(x$rho2_oof),
        rho2_iqr = stats::IQR(x$rho2_oof),
        fold_stable = mean(x$r2_oof > 0) >= 0.8 && median(x$rho2_oof) >= 0.25,
        stringsAsFactors = FALSE
    )
}))
rownames(locus_summary) <- NULL

cell_key <- interaction(locus_summary$population, locus_summary$region,
                        drop = TRUE, lex.order = TRUE)
cell_summary <- do.call(rbind, lapply(split(locus_summary, cell_key), function(x) {
    data.frame(
        population = x$population[[1L]],
        region = x$region[[1L]],
        selected_loci = nrow(x),
        stable_loci = sum(x$fold_stable),
        stable_locus_fraction = mean(x$fold_stable),
        diagnostic_pass = mean(x$fold_stable) >= 0.8,
        stringsAsFactors = FALSE
    )
}))
rownames(cell_summary) <- NULL

related <- audit[audit$sensitivity == "exclude_Br1105", , drop = FALSE]
relatedness_summary <- data.frame()
if (nrow(related)) {
    paired <- merge(
        baseline[, c("population_planned", "region_planned", "task_id_planned",
                     "seed_repeat", "rho2_oof", "r2_oof")],
        related[, c("population_planned", "region_planned", "task_id_planned",
                    "seed_repeat", "rho2_oof", "r2_oof")],
        by = c("population_planned", "region_planned", "task_id_planned",
               "seed_repeat"), suffixes = c("_full", "_exclude")
    )
    paired$rho2_difference <- paired$rho2_oof_exclude - paired$rho2_oof_full
    pair_key <- interaction(paired$population_planned, paired$region_planned,
                            drop = TRUE, lex.order = TRUE)
    relatedness_summary <- do.call(rbind, lapply(split(paired, pair_key), function(x) {
        data.frame(
            population = x$population_planned[[1L]],
            region = x$region_planned[[1L]],
            paired_runs = nrow(x),
            median_rho2_difference = median(x$rho2_difference),
            max_absolute_rho2_difference = max(abs(x$rho2_difference)),
            positive_r2_fraction_after_exclusion = mean(x$r2_oof_exclude > 0),
            stringsAsFactors = FALSE
        )
    }))
    rownames(relatedness_summary) <- NULL
}

dir.create(cli$output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(audit, file.path(cli$output_dir, "high-signal-audit-all-runs.tsv"))
write_tsv(locus_summary, file.path(cli$output_dir, "fold-stability-by-locus.tsv"))
write_tsv(cell_summary, file.path(cli$output_dir, "fold-stability-by-cell.tsv"))
write_tsv(relatedness_summary,
          file.path(cli$output_dir, "relatedness-sensitivity-by-cell.tsv"))
cat("High-signal diagnostic cells passing:", sum(cell_summary$diagnostic_pass),
    "of", nrow(cell_summary), "\n")
