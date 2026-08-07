#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

ROOT <- "/projects/b1213/users/kynon/projects/dna-methylation-heritability"
source(file.path(ROOT, "inputs/cell_proportions/_h/dnam_deconvolution_utils.R"))
cfg <- yaml::read_yaml(file.path(ROOT, "config/cell_deconvolution.yml"))
render <- function(template, region) gsub("\\{region\\}", region, template)
regions <- c("dlpfc", "hippocampus")
out_root <- file.path(ROOT, cfg$paths$output_root)

broad <- function(x, reference_key = NULL) {
  x[, broad_class := fcase(
    cell_type %in% c("Excitatory_neuron", "Inhibitory_neuron"), "Total_neuron",
    cell_type == "Astrocyte", "Astrocyte",
    cell_type == "Microglia", "Microglia",
    cell_type == "Oligodendrocyte", "Oligodendrocyte",
    cell_type == "OPC", "OPC",
    cell_type == "Endothelial", "Vascular_exploratory",
    default = NA_character_
  )]
  by_cols <- c("region", "sample_id", reference_key, "broad_class")
  x[!is.na(broad_class), .(proportion = sum(proportion)), by = by_cols]
}

rna_broad <- function(x, reg) {
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
    by = .(sample_id, broad_class)][, region := reg]
}

correlate <- function(matched, extra_by = character()) {
  by_cols <- c("region", extra_by, "broad_class")
  ans <- matched[, {
    ok <- is.finite(proportion_dnam) & is.finite(proportion_rna)
    if (sum(ok) >= 10L) {
      test <- suppressWarnings(cor.test(proportion_dnam[ok], proportion_rna[ok],
                                       method = "spearman", exact = FALSE))
      list(n = sum(ok), rho = unname(test$estimate), p_value = test$p.value)
    } else list(n = sum(ok), rho = NA_real_, p_value = NA_real_)
  }, by = by_cols]
  ans[, fdr := p.adjust(p_value, method = "BH")]
  ans
}

dnam <- rbindlist(lapply(regions, function(reg) {
  fread(file.path(ROOT, render(cfg$paths$dnam_wgbs_proportions_template, reg)))
}))
rna <- rbindlist(lapply(regions, function(reg) {
  rna_broad(fread(file.path(ROOT, render(cfg$paths$rna_proportions_template, reg))), reg)
}))
dnam_broad <- broad(copy(dnam))
matched <- merge(dnam_broad, rna, by = c("region", "sample_id", "broad_class"),
                 suffixes = c("_dnam", "_rna"))
cors <- correlate(matched)
cors[broad_class == "Total_neuron", neuron_fdr := p.adjust(p_value, method = "BH")]

components <- rbindlist(lapply(regions, function(reg) {
  work <- file.path(ROOT, render(cfg$paths$dnam_wgbs_work_template, reg))
  fread(file.path(work, "dnam-scmd-component-proportions.tsv.gz"))
}), fill = TRUE)
components[, reference_name := sub("_.*$", "", component)]
ref_est <- components[, .(proportion = median(proportion, na.rm = TRUE)),
                      by = .(region, sample_id, reference_name, cell_type)]
ref_est[, proportion := proportion / sum(proportion),
        by = .(region, sample_id, reference_name)]
ref_broad <- broad(ref_est, "reference_name")
ref_matched <- merge(ref_broad, rna,
  by = c("region", "sample_id", "broad_class"), suffixes = c("_dnam", "_rna"))
ref_cors <- correlate(ref_matched, "reference_name")

fwrite(matched, file.path(out_root, "dnam-scmd-wgbs-reference-rna-matched.tsv.gz"), sep = "\t")
fwrite(cors, file.path(out_root, "dnam-scmd-wgbs-reference-rna-concordance.tsv"), sep = "\t")
fwrite(ref_cors, file.path(out_root, "dnam-scmd-wgbs-reference-sensitivity.tsv"), sep = "\t")
primary_cor <- fread(file.path(out_root, "dnam-scmd-rna-concordance.tsv"))
platform_comparison <- rbindlist(list(
  primary_cor[region %in% regions & broad_class == "Total_neuron",
    .(region, platform = "850K", n, rho, p_value, neuron_fdr)],
  cors[broad_class == "Total_neuron",
    .(region, platform = "coordinate_WGBS", n, rho, p_value, neuron_fdr)]
))
fwrite(platform_comparison, file.path(out_root, "dnam-scmd-platform-neuronal-comparison.tsv"), sep = "\t")
print(cors[broad_class == "Total_neuron"])
