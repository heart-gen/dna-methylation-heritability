## Plot Fisher's enrichment results for AA non-heritable VMR subgroups.
##
## Subgroup labels:
##   A — AA non-heritable, EA heritable      (LD-misclassification candidate)
##   B — AA non-heritable, EA non-heritable  (jointly non-heritable; existing analysis)
##   C — AA non-heritable, EA low-pred/absent (AA-specific)
##
## Key plot: education enrichment should be in B (and possibly C), NOT A.
## If A lacks education enrichment, LD-confounding concern is directly answered.
##
## Outputs:
##   subgroup_tileplot_logit.pdf
##   subgroup_education_dotplot.pdf
##   subgroup_counts_barplot.pdf

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ggplot2)
  library(scales)
})

## Main
OUT_DIR <- here(
  "environmental-analysis/all_individuals/tissue_compare/non_heritable_subgroup/_m"
)

# Load data

enrich <- fread(file.path(OUT_DIR, "subgroup_enrichment_analysis.tsv"))
counts <- fread(file.path(OUT_DIR, "subgroup_counts.tsv"))

subgroup_labels <- c(
  A = "A  (AA nh, EA heritable)",
  B = "B  (jointly non-heritable)",
  C = "C  (AA-specific)"
)

enrich[, sig        := FDR < 0.05 & !is.na(FDR)]
enrich[, Subgroup   := factor(Subgroup, levels = c("A", "B", "C"))]
enrich[, Tissue     := factor(Tissue, levels = c("Caudate", "DLPFC", "Hippocampus"))]

enrich_logit <- enrich[Test == "Logit"]

# Readable env labels
env_labels <- c(
  smoking          = "Smoking",
  nicotine         = "Nicotine",
  ethanol          = "Alcohol",
  antipsychotics   = "Antipsychotics",
  cocaine          = "Cocaine",
  codeine          = "Codeine",
  morphine         = "Morphine",
  amphetamines     = "Amphetamines",
  hx_sexual_abuse  = "Sexual abuse",
  hx_physical_abuse = "Physical abuse",
  less_than_hs     = "Education < HS",
  more_than_hs     = "Education > HS",
  single           = "Marital: single",
  previously_married = "Marital: prev. married"
)
enrich_logit[, Env_label := ifelse(Env %in% names(env_labels),
                                    env_labels[Env], Env)]

# Order env by maximum absolute effect across subgroups
env_order <- enrich_logit[
  , .(max_effect = max(abs(log2(OR + 1e-3)), na.rm = TRUE)), by = Env_label
][order(-max_effect), Env_label]
enrich_logit[, Env_label := factor(Env_label, levels = rev(env_order))]


# Tile plot

p_tile <- ggplot(
  enrich_logit,
  aes(x = Subgroup, y = Env_label, fill = log2(OR))
) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(sig, "*", "")), size = 3.5, vjust = 0.8) +
  facet_wrap(~Tissue, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0,
    limits = log2(c(0.05, 20)), oob = scales::squish,
    name = "log\u2082(OR)"
  ) +
  scale_x_discrete(labels = c(
    A = "A\n(AA nh,\nEA heritable)",
    B = "B\n(jointly\nnon-heritable)",
    C = "C\n(AA-\nspecific)"
  )) +
  labs(
    title = "Environmental enrichment within AA non-heritable VMR subgroups",
    subtitle = "Logit test  |  * FDR < 0.05  |  background = all AA non-heritable VMRs",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text      = element_text(face = "bold"),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(OUT_DIR, "subgroup_tileplot_logit.pdf"),
       p_tile, width = 11, height = 7)


# Education dot plot

edu <- enrich_logit[Env %in% c("less_than_hs", "more_than_hs")]
edu[, Edu_label := ifelse(Env == "less_than_hs", "< High school", "> High school")]

p_edu <- ggplot(edu,
                aes(x = Subgroup, y = log2(OR),
                    colour = sig, shape = Edu_label,
                    group = Edu_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  facet_wrap(~Tissue, nrow = 1) +
  scale_colour_manual(
    values = c("TRUE" = "#D6604D", "FALSE" = "grey65"),
    labels = c("TRUE" = "FDR < 0.05", "FALSE" = "NS"),
    name   = NULL
  ) +
  scale_shape_manual(values = c(16, 17), name = "Education level") +
  scale_x_discrete(labels = c(
    A = "A\n(AA nh, EA her.)",
    B = "B\n(jointly nh)",
    C = "C\n(AA-specific)"
  )) +
  labs(
    title    = "Education enrichment across AA non-heritable VMR subgroups",
    subtitle = "Logit test  |  key test: education should be enriched in B, absent in A",
    x = "Subgroup", y = "log\u2082(Odds Ratio)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text        = element_text(face = "bold"),
    panel.grid.minor  = element_blank(),
    legend.position   = "right"
  )

ggsave(file.path(OUT_DIR, "subgroup_education_dotplot.pdf"),
       p_edu, width = 9, height = 5)


# Subgroup count bar

counts[, subgroup := factor(subgroup, levels = c("B", "A", "C"))]
counts[, tissue   := factor(tissue, levels = c("Caudate", "DLPFC", "Hippocampus"))]

fill_cols <- c(
  B = "#4DAF4A",  # jointly non-heritable — green
  A = "#E41A1C",  # AA nh / EA heritable  — red
  C = "#377EB8"   # AA-specific           — blue
)
fill_labels <- c(
  B = "B: jointly non-heritable",
  A = "A: AA nh, EA heritable",
  C = "C: AA-specific"
)

p_counts <- ggplot(counts,
                   aes(x = tissue, y = n_vmrs,
                       fill = subgroup)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = fill_cols, labels = fill_labels, name = "Subgroup") +
  labs(
    title = "AA non-heritable VMR subgroup sizes by brain region",
    x = NULL, y = "Number of VMRs"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position      = "right",
    panel.grid.major.x   = element_blank()
  )

ggsave(file.path(OUT_DIR, "subgroup_counts_barplot.pdf"),
       p_counts, width = 7, height = 5)

message("All plots saved to: ", OUT_DIR)
sessionInfo()
