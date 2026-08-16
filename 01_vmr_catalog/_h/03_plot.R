#### 01_vmr_catalog / 03_plot: catalog diagnostics ####
##
## Diagnostic figures only. Manuscript figures are assembled exclusively in
## 09_integrated_manuscript_outputs from accepted run IDs (AGENTS.md 7.9).
##
## Usage:
##   Rscript 03_plot.R --cohort AA --region caudate --run-id ID

source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

opts <- parse_v2_args(require = c("cohort", "region", "run_id"))
cohort <- opts$cohort; region <- opts$region

module_root <- file.path(V2_ROOT, "01_vmr_catalog")
run_dir <- file.path(module_root, "_m", "runs", opts$run_id)
vmr_dir <- file.path(run_dir, "vmr")
fig_dir <- file.path(run_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

vmr <- fread(file.path(vmr_dir, "vmr_catalog.tsv"))
cutoffs <- fread(file.path(vmr_dir, "sd_cutoffs.tsv"))
membership <- fread(file.path(vmr_dir, "cpg_vmr_membership.tsv"))

tag <- paste0(cohort, " / ", region)
theme_set(theme_bw(base_size = 11))

chr_levels <- paste0("chr", chrom_order(include_sex = FALSE))
ord <- function(x) factor(x, levels = chr_levels)

## VMR count per chromosome
ggsave(file.path(fig_dir, "vmr_count_by_chrom.pdf"),
       ggplot(cutoffs, aes(ord(chr), n_vmrs)) +
           geom_col(fill = "grey30") +
           labs(x = NULL, y = "VMRs", title = paste("VMRs per chromosome —", tag)) +
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5)),
       width = 8, height = 4)

## Per-chromosome SD cutoff. Worth watching: the audit found the V1 misalignment
## inflated this cutoff by ~9.2%, and colSds() uses n-1 on residuals from a
## 6-parameter fit, so cutoffs are not directly comparable across regions.
ggsave(file.path(fig_dir, "sd_cutoff_by_chrom.pdf"),
       ggplot(cutoffs, aes(ord(chr), sd_cutoff)) +
           geom_point(size = 2) +
           labs(x = NULL, y = "residual SD cutoff (99th pct)",
                title = paste("SD cutoff per chromosome —", tag)) +
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5)),
       width = 8, height = 4)

## VMR width and CpG-count distributions
vmr[, width := end - start + 1L]
ggsave(file.path(fig_dir, "vmr_width.pdf"),
       ggplot(vmr, aes(width)) +
           geom_histogram(bins = 60, fill = "grey30") +
           scale_x_log10() +
           labs(x = "VMR width (bp, log10)", y = "count",
                title = paste("VMR width —", tag)),
       width = 6, height = 4)

ggsave(file.path(fig_dir, "vmr_cpg_count.pdf"),
       ggplot(vmr, aes(n)) +
           geom_histogram(bins = 60, fill = "grey30") +
           scale_x_log10() +
           labs(x = "CpGs per VMR (log10)", y = "count",
                title = paste("CpGs per VMR —", tag)),
       width = 6, height = 4)

## Residual SD of constituent CpGs
ggsave(file.path(fig_dir, "vmr_cpg_residual_sd.pdf"),
       ggplot(membership, aes(cpg_residual_sd)) +
           geom_histogram(bins = 60, fill = "grey30") +
           labs(x = "CpG residual SD", y = "count",
                title = paste("Residual SD of VMR CpGs —", tag)),
       width = 6, height = 4)

## Source data for every panel, so a figure can be traced to numbers.
write_atomic(cutoffs, file.path(fig_dir, "source_data_by_chrom.tsv"))
write_atomic(vmr[, .(vmr_id, chr, start, end, n, width)],
             file.path(fig_dir, "source_data_vmr.tsv"))

message("[done] figures written to ", fig_dir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
