#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(yaml)
})

ROOT <- "/projects/b1213/users/kynon/projects/dna-methylation-heritability"
cfg <- yaml::read_yaml(file.path(ROOT, "config/cell_deconvolution.yml"))
render <- function(template, region) gsub("\\{region\\}", region, template)
out_root <- file.path(ROOT, cfg$paths$output_root)

props <- rbindlist(lapply(cfg$regions, function(region) {
  fread(file.path(ROOT, render(cfg$paths$dnam_proportions_template, region)))
}))
matched <- fread(file.path(out_root, "dnam-scmd-rna-matched.tsv.gz"))
cors <- fread(file.path(out_root, "dnam-scmd-rna-concordance.tsv"))
sqc <- rbindlist(lapply(cfg$regions, function(region) {
  fread(file.path(ROOT, render(cfg$paths$dnam_work_template, region), "sample_marker_qc.tsv"))
}))

region_labels <- c(caudate = "Caudate", dlpfc = "DLPFC", hippocampus = "Hippocampus")
cell_order <- c("Excitatory_neuron", "Inhibitory_neuron", "Astrocyte", "Oligodendrocyte",
                "OPC", "Microglia", "Endothelial")
props[, cell_type := factor(cell_type, levels = cell_order)]

p1 <- ggplot(props, aes(cell_type, proportion, fill = cell_type)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25, fill = "white") +
  facet_wrap(~factor(region, levels = names(region_labels), labels = region_labels), nrow = 1) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "DNAm-estimated proportion", title = "A  DNAm cell-composition estimates") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

neur <- matched[broad_class == "Total_neuron"]
p2 <- ggplot(neur, aes(proportion_rna, proportion_dnam)) +
  geom_point(alpha = 0.55, size = 1.1, color = "#2C7FB8") +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, color = "#D95F0E") +
  facet_wrap(~factor(region, levels = names(region_labels), labels = region_labels), scales = "free") +
  labs(x = "RNA MuSiC total-neuronal proportion", y = "DNAm total-neuronal proportion",
       title = "B  Cross-modality neuronal concordance") +
  theme_bw(base_size = 9)

p3 <- ggplot(cors, aes(factor(region, levels = names(region_labels), labels = region_labels),
                       broad_class, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(is.finite(rho), sprintf("%.2f", rho), "NA")), size = 2.6) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0,
                       limits = c(-1, 1), name = "Spearman rho") +
  labs(x = NULL, y = NULL, title = "C  DNAm–RNA cell-class correlations") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

p4 <- ggplot(sqc, aes(marker_coverage_fraction, fill = region)) +
  geom_histogram(binwidth = 0.01, boundary = 0, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = as.numeric(cfg$filtering$min_sample_marker_fraction), linetype = 2) +
  facet_wrap(~factor(region, levels = names(region_labels), labels = region_labels), nrow = 1) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Fraction of retained markers with coverage >=5", y = "Samples",
       title = "D  Sample marker coverage") +
  theme_bw(base_size = 9) + theme(legend.position = "none")

fig <- (p1 / p2 / (p3 | p4)) + plot_layout(heights = c(1.1, 1, 0.9))
ggsave(file.path(out_root, "dnam_scmd_qc_figure.pdf"), fig, width = 12, height = 11,
       device = cairo_pdf)
ggsave(file.path(out_root, "dnam_scmd_qc_figure.png"), fig, width = 12, height = 11,
       dpi = 300)
