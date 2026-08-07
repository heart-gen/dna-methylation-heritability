#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

ROOT <- "/projects/b1213/users/kynon/projects/dna-methylation-heritability"
source(file.path(ROOT, "inputs/cell_proportions/_h/dnam_deconvolution_utils.R"))
cfg <- yaml::read_yaml(file.path(ROOT, "config/cell_deconvolution.yml"))
render <- function(template, region) gsub("\\{region\\}", region, template)
out_root <- file.path(ROOT, cfg$paths$output_root)

to_broad_dnam <- function(x) {
  x[, broad_class := fcase(
    cell_type %in% c("Excitatory_neuron", "Inhibitory_neuron"), "Total_neuron",
    cell_type == "Astrocyte", "Astrocyte",
    cell_type == "Microglia", "Microglia",
    cell_type == "Oligodendrocyte", "Oligodendrocyte",
    cell_type == "OPC", "OPC",
    cell_type == "Endothelial", "Vascular_exploratory",
    default = NA_character_
  )]
  x[!is.na(broad_class), .(proportion = sum(proportion)),
    by = .(region, sample_id, broad_class)]
}

to_broad_rna <- function(x, region) {
  x[, broad_class := fcase(
    cell_type %in% c("Excit", "Inhib", "D1-SPN", "D2-SPN"), "Total_neuron",
    cell_type == "Astro", "Astrocyte",
    cell_type == "Micro", "Microglia",
    cell_type == "Oligo", "Oligodendrocyte",
    cell_type == "OPC", "OPC",
    cell_type == "Mural", "Vascular_exploratory",
    default = NA_character_
  )]
  x[!is.na(broad_class), .(proportion = sum(proportion)),
    by = .(sample_id, broad_class)][, region := region]
}

all_dnam <- list()
all_rna <- list()
sample_qc <- list()
fraction_qc <- list()
component_props <- list()
for (region in cfg$regions) {
  dfile <- file.path(ROOT, render(cfg$paths$dnam_proportions_template, region))
  rfile <- file.path(ROOT, render(cfg$paths$rna_proportions_template, region))
  if (!file.exists(dfile)) stop("Missing DNAm proportions: ", dfile)
  if (!file.exists(rfile)) stop("Missing RNA proportions: ", rfile)
  all_dnam[[region]] <- to_broad_dnam(fread(dfile))
  all_rna[[region]] <- to_broad_rna(fread(rfile), region)
  work <- file.path(ROOT, render(cfg$paths$dnam_work_template, region))
  sample_qc[[region]] <- fread(file.path(work, "sample_marker_qc.tsv"))
  fraction_qc[[region]] <- fread(file.path(work, "deconvolution_qc.tsv"))
  component_props[[region]] <- fread(file.path(work, "dnam-scmd-component-proportions.tsv.gz"))
}

dnam <- rbindlist(all_dnam)
rna <- rbindlist(all_rna)
matched <- merge(dnam, rna, by = c("region", "sample_id", "broad_class"),
                 suffixes = c("_dnam", "_rna"))

cor_rows <- matched[, {
  ok <- is.finite(proportion_dnam) & is.finite(proportion_rna)
  if (sum(ok) >= 10L) {
    test <- suppressWarnings(cor.test(proportion_dnam[ok], proportion_rna[ok],
                                     method = "spearman", exact = FALSE))
    list(n = sum(ok), rho = unname(test$estimate), p_value = test$p.value)
  } else list(n = sum(ok), rho = NA_real_, p_value = NA_real_)
}, by = .(region, broad_class)]
cor_rows[, fdr_all := p.adjust(p_value, method = "BH")]
cor_rows[broad_class == "Total_neuron", neuron_fdr := p.adjust(p_value, method = "BH")]

# Lee/Tian sensitivity: take the median across component algorithms within each
# reference, re-close fractions, and repeat broad cross-modality correlations.
component <- rbindlist(component_props, fill = TRUE)
component[, reference_name := sub("_.*$", "", component)]
reference_estimates <- component[, .(proportion = median(proportion, na.rm = TRUE)),
                                 by = .(region, sample_id, reference_name, cell_type)]
reference_estimates[, proportion := proportion / sum(proportion),
                    by = .(region, sample_id, reference_name)]
# Re-aggregate directly to retain the reference key (to_broad_dnam intentionally
# emits only the primary output keys).
reference_estimates[, broad_class := fcase(
  cell_type %in% c("Excitatory_neuron", "Inhibitory_neuron"), "Total_neuron",
  cell_type == "Astrocyte", "Astrocyte",
  cell_type == "Microglia", "Microglia",
  cell_type == "Oligodendrocyte", "Oligodendrocyte",
  cell_type == "OPC", "OPC",
  cell_type == "Endothelial", "Vascular_exploratory",
  default = NA_character_
)]
reference_broad <- reference_estimates[!is.na(broad_class),
  .(proportion = sum(proportion)),
  by = .(region, sample_id, reference_name, broad_class)]
reference_matched <- merge(reference_broad, rna,
  by = c("region", "sample_id", "broad_class"), suffixes = c("_dnam", "_rna"))
reference_cor <- reference_matched[, {
  ok <- is.finite(proportion_dnam) & is.finite(proportion_rna)
  if (sum(ok) >= 10L) {
    test <- suppressWarnings(cor.test(proportion_dnam[ok], proportion_rna[ok],
                                     method = "spearman", exact = FALSE))
    list(n = sum(ok), rho = unname(test$estimate), p_value = test$p.value)
  } else list(n = sum(ok), rho = NA_real_, p_value = NA_real_)
}, by = .(region, reference_name, broad_class)]
reference_cor[, fdr := p.adjust(p_value, method = "BH")]

phen <- fread(file.path(ROOT, cfg$paths$phenotype_table))
phen[, region := normalize_region(region)]
phen <- phen[, .(sample_id = as.character(brnum), region, pmi, ph)]
sqc <- rbindlist(sample_qc, fill = TRUE)
tech_base <- merge(dnam, sqc[, .(sample_id, region, mean_marker_coverage)],
                   by = c("sample_id", "region"), all.x = TRUE)
tech_base <- merge(tech_base, phen, by = c("sample_id", "region"), all.x = TRUE)
technical <- rbindlist(lapply(c("mean_marker_coverage", "pmi", "ph"), function(variable) {
  tech_base[, {
    y <- get(variable)
    ok <- is.finite(proportion) & is.finite(y)
    if (sum(ok) >= 10L) {
      test <- suppressWarnings(cor.test(proportion[ok], y[ok], method = "spearman", exact = FALSE))
      list(variable = variable, n = sum(ok), rho = unname(test$estimate), p_value = test$p.value)
    } else list(variable = variable, n = sum(ok), rho = NA_real_, p_value = NA_real_)
  }, by = .(region, broad_class)]
}))
technical[, fdr := p.adjust(p_value, method = "BH")]
technical_availability <- data.table(
  covariate = c("mean_marker_coverage", "pmi", "ph", "processing_batch"),
  available = c(TRUE, "pmi" %in% names(phen), "ph" %in% names(phen), FALSE),
  source = c("WGBS marker extraction", "phenotype_data.tsv", "phenotype_data.tsv", NA_character_),
  note = c(
    "tested",
    "tested",
    "tested",
    "No processing-batch field was present in phenotype_data.tsv or the HDF5 BSseq colData"
  )
)

fqc <- rbindlist(fraction_qc, fill = TRUE)
validation <- rbindlist(lapply(cfg$regions, function(reg) {
  nrow_expected <- uniqueN(phen[region == reg]$sample_id)
  n_est <- uniqueN(fqc[region == reg]$sample_id)
  neuron <- cor_rows[region == reg & broad_class == "Total_neuron"]
  bound_ok <- nrow(fqc[region == reg]) > 0L && all(fqc[region == reg]$bounded)
  sum_ok <- nrow(fqc[region == reg]) > 0L && all(fqc[region == reg]$sum_to_one)
  sample_fraction <- n_est / nrow_expected
  neuron_ok <- nrow(neuron) == 1L && is.finite(neuron$rho) && is.finite(neuron$neuron_fdr) &&
    neuron$rho >= as.numeric(cfg$validation$min_neuronal_spearman_rho) &&
    neuron$neuron_fdr <= as.numeric(cfg$validation$max_neuronal_fdr)
  data.table(
    region = reg,
    n_expected_samples = nrow_expected,
    n_estimated_samples = n_est,
    estimated_sample_fraction = sample_fraction,
    all_fractions_bounded = bound_ok,
    all_fractions_sum_to_one = sum_ok,
    neuronal_rho = if (nrow(neuron)) neuron$rho else NA_real_,
    neuronal_fdr = if (nrow(neuron)) neuron$neuron_fdr else NA_real_,
    integration_pass = bound_ok && sum_ok &&
      sample_fraction >= 1 - as.numeric(cfg$filtering$max_failed_sample_fraction_for_integration) &&
      neuron_ok
  )
}))

fwrite(matched, file.path(out_root, "dnam-scmd-rna-matched.tsv.gz"), sep = "\t")
fwrite(cor_rows, file.path(out_root, "dnam-scmd-rna-concordance.tsv"), sep = "\t")
fwrite(reference_estimates, file.path(out_root, "dnam-scmd-reference-estimates.tsv.gz"), sep = "\t")
fwrite(reference_cor, file.path(out_root, "dnam-scmd-reference-sensitivity.tsv"), sep = "\t")
fwrite(technical, file.path(out_root, "dnam-scmd-technical-correlations.tsv"), sep = "\t")
fwrite(technical_availability, file.path(out_root, "dnam-scmd-technical-covariate-availability.tsv"), sep = "\t")
fwrite(validation, file.path(out_root, "dnam-scmd-validation-summary.tsv"), sep = "\t")
for (reg in cfg$regions) {
  fwrite(validation[region == reg],
         file.path(out_root, paste0("dnam-scmd-validation-", reg, ".tsv")), sep = "\t")
}
capture.output(sessionInfo(), file = file.path(out_root, "dnam_scmd_validation_session_info.txt"))
print(validation)
