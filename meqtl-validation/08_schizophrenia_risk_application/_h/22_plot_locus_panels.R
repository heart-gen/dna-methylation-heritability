#!/usr/bin/env Rscript
# Manuscript-quality SCZ locus panels.
# Inputs from 21_prepare_locus_panel_data.py under _m/locus_panels/<rsID>/

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
})

root <- "/projects/b1213/users/kynon/projects/dna-methylation-heritability"
base <- file.path(root, "meqtl-validation/08_schizophrenia_risk_application/_m/locus_panels")
fig_dir <- file.path(base, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

region_cols <- c(
  caudate = "#1B4F72",
  dlpfc = "#B9770E",
  hippocampus = "#196F3D"
)

theme_locus <- function(base_size = 9) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1, color = "black"),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.margin = margin(4, 6, 4, 6),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size)
    )
}

plot_one_locus <- function(snp) {
  d <- file.path(base, snp)
  meta <- fread(file.path(d, "locus_meta.tsv"))
  gwas <- fread(file.path(d, "gwas_regional_hg19.tsv.gz"))
  meq <- fread(file.path(d, "meqtl_caudate_cpgs.tsv.gz"))
  forest <- fread(file.path(d, "cross_region_forest.tsv.gz"))
  vmr <- fread(file.path(d, "vmr_predictability.tsv"))
  tx <- tryCatch(fread(file.path(d, "tx_links_fdr.tsv")), error = function(e) data.table())
  gtex <- tryCatch(fread(file.path(d, "gtex_level3_hits.tsv")), error = function(e) data.table())

  pos0 <- meta$pos_hg19[1]
  gwas[, neglog10p := -log10(pmax(pval, 1e-300))]
  gwas[, is_index := snp == meta$index_snp[1] | abs(pos_hg19 - pos0) < 2]

  # A: GWAS
  pA <- ggplot(gwas, aes(x = pos_hg19 / 1e6, y = neglog10p)) +
    geom_point(aes(color = is_index), size = 0.7, alpha = 0.75, stroke = 0) +
    scale_color_manual(values = c("FALSE" = "grey55", "TRUE" = "#8B0000"), guide = "none") +
    geom_point(
      data = gwas[is_index == TRUE],
      aes(x = pos_hg19 / 1e6, y = neglog10p),
      color = "#8B0000", size = 2.2
    ) +
    labs(x = paste0("chr", meta$chrom_hg19[1], " position (Mb, hg19)"), y = expression(-log[10](P)[GWAS])) +
    theme_locus()

  # B: caudate index-SNP meQTL along CpGs
  if (nrow(meq) > 0) {
    meq[, neglog10p := -log10(pmax(pval_nominal, 1e-300))]
    meq[, sig := significant_fdr == TRUE | significant_fdr == "True" | significant_fdr == "TRUE"]
    # VMR intervals in hg38 CpG coords — plot vs cpg_pos
    vmr_rects <- vmr[!is.na(start) & !is.na(end)]
    pB <- ggplot(meq, aes(x = cpg_pos / 1e6, y = neglog10p)) +
      {
        if (nrow(vmr_rects) > 0) {
          geom_rect(
            data = vmr_rects,
            aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf, fill = is_best),
            inherit.aes = FALSE, alpha = 0.12, color = NA
          )
        }
      } +
      geom_point(aes(color = sig), size = 1.3, alpha = 0.9) +
      scale_color_manual(
        values = c("FALSE" = "grey50", "TRUE" = "#1B4F72", "False" = "grey50", "True" = "#1B4F72"),
        labels = c("FALSE" = "NS", "TRUE" = "FDR<0.05", "False" = "NS", "True" = "FDR<0.05"),
        name = NULL
      ) +
      scale_fill_manual(values = c("TRUE" = "#1B4F72", "FALSE" = "#85C1E9", "True" = "#1B4F72", "False" = "#85C1E9"), guide = "none") +
      labs(x = paste0(meta$chrom_hg38[1], " CpG position (Mb, hg38)"), y = expression(-log[10](P)[meQTL])) +
      theme_locus() +
      theme(legend.position = c(0.85, 0.85), legend.key.size = unit(0.3, "cm"))
  } else {
    pB <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "No meQTL tests")
  }

  # C: cross-region forest (top CpGs)
  if (nrow(forest) > 0) {
    forest[, region := factor(region, levels = c("caudate", "dlpfc", "hippocampus"))]
    forest[, phenotype_id := factor(phenotype_id, levels = unique(phenotype_id))]
    forest[, lo := beta - 1.96 * se]
    forest[, hi := beta + 1.96 * se]
    pC <- ggplot(forest, aes(x = beta, y = phenotype_id, color = region)) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.3) +
      geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, linewidth = 0.6, orientation = "y", position = position_dodge(width = 0.6)) +
      geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
      scale_color_manual(values = region_cols, name = NULL) +
      labs(x = "meQTL beta (risk allele)", y = NULL) +
      theme_locus() +
      theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
  } else {
    pC <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "No forest pairs")
  }

  # D: VMR predictability + compact multi-omic tags
  if (nrow(vmr) > 0) {
    vmr[, vmr_lab := paste0("VMR ", vmr_id)]
    vmr[, vmr_lab := factor(vmr_lab, levels = vmr_lab[order(local_predictability)])]
    tx_bits <- character(0)
    if (isTRUE(as.logical(meta$tx_expression[1]))) tx_bits <- c(tx_bits, "expression")
    if (isTRUE(as.logical(meta$tx_psi[1]))) tx_bits <- c(tx_bits, "PSI")
    ann_lines <- character(0)
    if (length(tx_bits) > 0) {
      ann_lines <- c(ann_lines, paste0("VMR-transcript FDR: ", paste(tx_bits, collapse = " + ")))
    } else {
      ann_lines <- c(ann_lines, "VMR-transcript FDR: none")
    }
    if (nrow(tx) > 0 && "gene_symbol" %in% names(tx)) {
      tx_genes <- unique(sprintf("%s (%s)", tx$gene_symbol, tx$modality))
      ann_lines <- c(ann_lines, paste(head(tx_genes, 3), collapse = ", "))
    }
    if (isTRUE(as.logical(meta$level3_pass[1]))) {
      gtex_lab <- if (nrow(gtex) > 0 && "gene_symbol" %in% names(gtex)) {
        paste(head(unique(as.character(gtex$gene_symbol[!is.na(gtex$gene_symbol) & gtex$gene_symbol != ""])), 4), collapse = ", ")
      } else {
        "yes"
      }
      ann_lines <- c(ann_lines, paste0("Level 3 (GTEx eQTL): ", gtex_lab))
    } else {
      ann_lines <- c(ann_lines, "Level 3 (GTEx eQTL): no")
    }
    ann <- paste(ann_lines, collapse = "\n")
    xmax <- max(vmr$local_predictability, na.rm = TRUE)
    vmr[, fill_lab := ifelse(
      as.logical(is_best) %in% TRUE,
      "Index meQTL VMR",
      "Other linked VMR"
    )]
    pD <- ggplot(vmr, aes(x = local_predictability, y = vmr_lab, fill = fill_lab)) +
      geom_col(width = 0.65) +
      geom_text(
        aes(label = sprintf("%.2f", local_predictability)),
        hjust = -0.15, size = 2.6, color = "grey20"
      ) +
      scale_fill_manual(
        values = c("Index meQTL VMR" = "#1B4F72", "Other linked VMR" = "#AED6F1"),
        name = NULL
      ) +
      scale_x_continuous(limits = c(0, xmax * 1.65), expand = expansion(mult = c(0, 0.02))) +
      labs(x = "Local genetic predictability", y = NULL) +
      annotate(
        "text", x = xmax * 1.62, y = Inf, label = ann,
        hjust = 1, vjust = 1.15, size = 2.5, color = "grey25", lineheight = 1.05
      ) +
      theme_locus() +
      theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
  } else {
    pD <- ggplot() + theme_void()
  }

  tag_theme <- theme(plot.tag = element_text(face = "bold", size = 11))
  fig <- (pA + labs(tag = "A") + tag_theme) +
    (pB + labs(tag = "B") + tag_theme) +
    (pC + labs(tag = "C") + tag_theme) +
    (pD + labs(tag = "D") + tag_theme) +
    plot_layout(ncol = 2, heights = c(1, 1.15))

  outfile <- file.path(fig_dir, paste0(snp, "_locus_panel.pdf"))
  ggsave(outfile, fig, width = 7.2, height = 6.2, units = "in", useDingbats = FALSE)
  ggsave(sub("\\.pdf$", ".png", outfile), fig, width = 7.2, height = 6.2, units = "in", dpi = 300)
  message("Wrote ", outfile)
  invisible(outfile)
}

manifest <- fread(file.path(base, "locus_panel_manifest.tsv"))
outs <- vapply(manifest$index_snp, plot_one_locus, character(1))

# Combined hero strip for main text (top 2 heroes + strongest predictability)
heroes <- manifest[hero == TRUE | rank <= 3][order(rank)]
if (nrow(heroes) >= 1) {
  # Also write an index README for manuscript
  sink(file.path(fig_dir, "LOCUS_PANEL_README.md"))
  cat("# SCZ locus panels\n\n")
  cat("Generated by `22_plot_locus_panels.R` from tidy tables in `../<rsID>/`.\n\n")
  cat("| Rank | Index SNP | PDF | Hero |\n|---:|---|---|---|\n")
  for (i in seq_len(nrow(manifest))) {
    cat(sprintf(
      "| %s | %s | `%s_locus_panel.pdf` | %s |\n",
      manifest$rank[i], manifest$index_snp[i], manifest$index_snp[i],
      ifelse(isTRUE(manifest$hero[i]), "yes", "no")
    ))
  }
  cat("\nPanel key: **A** PGC3 GWAS regional association (hg19); **B** caudate index-SNP CpG meQTL (hg38); **C** cross-region meQTL forest; **D** VMR local genetic predictability + TX/Level3 annotation.\n")
  cat("\nRecommended main-text heroes: rs8048039, rs13331198. Others supplemental.\n")
  sink()
}

message("Done. Figures in ", fig_dir)
