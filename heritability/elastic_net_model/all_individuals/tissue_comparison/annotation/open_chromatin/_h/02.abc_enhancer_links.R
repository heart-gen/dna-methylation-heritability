#### ABC Enhancer-Promoter Links for Intergenic VMRs ####
##
## Connects intergenic VMRs to putative target genes using BrainScope ABC
## enhancer-promoter links. Compares target gene sets between heritable and
## non-heritable intergenic VMRs, and runs gene set enrichment (GO/KEGG)
## using clusterProfiler.
##
## Input:  open_chromatin/_m/intergenic_vmr_atac_overlap.tsv
##         inputs/brainscope/_m/High_ABC_results_CPM1only_Concise_RK_7_8_21.csv
##
## Note: ABC data derived from PFC; applied to all 3 tissues with caveat.
##       clusterProfiler requires the 'rnaseq' conda environment.
##
## Run: conda run -p $ENV_PATH/rnaseq Rscript ../_h/02.abc_enhancer_links.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(GenomicRanges)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

## Configuration

BRAINSCOPE_DIR <- here::here("inputs", "brainscope", "_m")
IN_DIR <- here::here(
  "heritability", "elastic_net_model", "BA_only",
  "tissue_comparison", "annotation", "open_chromatin", "_m"
)
OUT_DIR <- IN_DIR

ABC_FILE <- file.path(BRAINSCOPE_DIR,
  "High_ABC_results_CPM1only_Concise_RK_7_8_21.csv")

N_QUINTILES <- 5

## Load intergenic VMR data (from script 01)

cat("Loading intergenic VMR ATAC overlap table...\n")
vmr_df <- fread(file.path(IN_DIR, "intergenic_vmr_atac_overlap.tsv"))
cat(sprintf("  Total intergenic VMRs: %d (%d tissues)\n",
  nrow(vmr_df), length(unique(vmr_df$tissue))))

# Pooled (unique VMRs across tissues by coordinates)
vmr_pooled <- vmr_df |>
  distinct(seqnames, start, end, h2_category, h2_unscaled, .keep_all = TRUE)
cat(sprintf("  Unique intergenic VMRs (pooled): %d\n", nrow(vmr_pooled)))

## Load ABC links

cat("\nLoading ABC enhancer-promoter links...\n")
abc <- fread(ABC_FILE)
cat(sprintf("  ABC links: %d rows\n", nrow(abc)))
cat(sprintf("  Columns: %s\n", paste(colnames(abc), collapse = ", ")))

# Create GRanges for ABC enhancer peaks
abc_gr <- GRanges(
  seqnames = abc$PeakID_seqnames,
  ranges   = IRanges(start = abc$PeakID_start, end = abc$PeakID_end)
)
mcols(abc_gr)$TargetGene_name <- abc$TargetGene_name
mcols(abc_gr)$ABC.Score       <- abc$ABC.Score
mcols(abc_gr)$distance        <- abc$distance
mcols(abc_gr)$TargetGene      <- abc$TargetGene
mcols(abc_gr)$abc_row         <- seq_len(nrow(abc))

## Overlap VMRs with ABC enhancer peaks

cat("\nOverlapping VMRs with ABC enhancer peaks...\n")

overlap_abc <- function(df) {
  vmr_gr <- GRanges(
    seqnames = df$seqnames,
    ranges   = IRanges(start = df$start, end = df$end)
  )
  hits <- findOverlaps(vmr_gr, abc_gr, ignore.strand = TRUE)
  if (length(hits) == 0) return(tibble())

  vmr_idx <- queryHits(hits)
  abc_idx <- subjectHits(hits)

  tibble(
    seqnames       = df$seqnames[vmr_idx],
    start          = df$start[vmr_idx],
    end            = df$end[vmr_idx],
    h2_category    = df$h2_category[vmr_idx],
    h2_unscaled    = df$h2_unscaled[vmr_idx],
    tissue         = df$tissue[vmr_idx],
    TargetGene_name = abc$TargetGene_name[abc_idx],
    TargetGene      = abc$TargetGene[abc_idx],
    ABC.Score       = abc$ABC.Score[abc_idx],
    distance        = abc$distance[abc_idx]
  )
}

# Per-tissue links
cat("  Per-tissue overlap...\n")
abc_links_per_tissue <- lapply(unique(vmr_df$tissue), function(t) {
  overlap_abc(vmr_df[vmr_df$tissue == t, ])
})
abc_links_tissue <- bind_rows(abc_links_per_tissue)
cat(sprintf("  Per-tissue VMR-ABC links: %d\n", nrow(abc_links_tissue)))

# Pooled links
cat("  Pooled overlap...\n")
abc_links_pooled <- overlap_abc(vmr_pooled)
abc_links_pooled$tissue <- "Pooled"
cat(sprintf("  Pooled VMR-ABC links: %d\n", nrow(abc_links_pooled)))

# Combine
abc_links <- bind_rows(abc_links_tissue, abc_links_pooled)
fwrite(abc_links, file.path(OUT_DIR, "abc_vmr_gene_links.tsv"), sep = "\t")
cat("Saved: abc_vmr_gene_links.tsv\n")

## Extract target gene lists per comparison group

get_genes <- function(df, h2_cat) {
  df |>
    filter(h2_category == h2_cat) |>
    pull(TargetGene_name) |>
    unique() |>
    sort()
}

# Per tissue
gene_lists <- lapply(c(unique(vmr_df$tissue), "Pooled"), function(t) {
  df_t <- abc_links[abc_links$tissue == t, ]
  list(
    Heritable    = get_genes(df_t, "Heritable"),
    NonHeritable = get_genes(df_t, "Non-heritable")
  )
})
names(gene_lists) <- c(unique(vmr_df$tissue), "Pooled")

# Save gene lists
for (t in names(gene_lists)) {
  t_safe <- gsub("[^A-Za-z0-9]", "_", t)
  writeLines(gene_lists[[t]]$Heritable,
    file.path(OUT_DIR, paste0("target_genes_heritable_", t_safe, ".txt")))
  writeLines(gene_lists[[t]]$NonHeritable,
    file.path(OUT_DIR, paste0("target_genes_nonheritable_", t_safe, ".txt")))
}
cat("Saved target gene list files.\n")

# Top quintile vs Q1-Q4 gene lists (pooled)
cat("\nComputing Q5 vs Q1-Q4 target gene lists...\n")
vmr_pooled2 <- vmr_pooled |>
  mutate(
    h2_quintile_cut = cut(
      h2_unscaled,
      breaks = quantile(h2_unscaled, probs = seq(0, 1, 0.2), na.rm = TRUE),
      labels = paste0("Q", 1:5),
      include.lowest = TRUE
    )
  )
abc_q <- overlap_abc(vmr_pooled2 |> mutate(tissue = "Pooled"))
abc_q2 <- abc_q |>
  left_join(
    vmr_pooled2 |> dplyr::select(seqnames, start, end, h2_quintile_cut),
    by = c("seqnames", "start", "end")
  )

genes_q5    <- abc_q2 |> filter(h2_quintile_cut == "Q5") |> pull(TargetGene_name) |> unique() |> sort()
genes_q1q4  <- abc_q2 |> filter(h2_quintile_cut != "Q5") |> pull(TargetGene_name) |> unique() |> sort()
writeLines(genes_q5,   file.path(OUT_DIR, "target_genes_q5_pooled.txt"))
writeLines(genes_q1q4, file.path(OUT_DIR, "target_genes_q1q4_pooled.txt"))
cat(sprintf("  Q5 target genes: %d | Q1-Q4 target genes: %d\n",
  length(genes_q5), length(genes_q1q4)))

## ABC score comparison: heritable vs non-heritable (pooled)

abc_score_summary <- abc_links |>
  filter(tissue == "Pooled") |>
  group_by(h2_category) |>
  summarise(
    n_links    = n(),
    n_vmrs     = n_distinct(paste0(seqnames, ":", start, "-", end)),
    n_genes    = n_distinct(TargetGene_name),
    median_abc = median(ABC.Score, na.rm = TRUE),
    mean_abc   = mean(ABC.Score,   na.rm = TRUE),
    median_dist = median(distance, na.rm = TRUE),
    .groups    = "drop"
  )
fwrite(abc_score_summary, file.path(OUT_DIR, "abc_score_summary.tsv"), sep = "\t")
cat("Saved: abc_score_summary.tsv\n")
print(abc_score_summary)

# Wilcoxon test: ABC score by h2 category
wt <- wilcox.test(
  ABC.Score ~ h2_category,
  data = abc_links |> filter(tissue == "Pooled",
                              h2_category %in% c("Heritable", "Non-heritable"))
)
cat(sprintf("\nWilcoxon test (ABC score, Heritable vs Non-heritable): W=%.0f, p=%.3g\n",
  wt$statistic, wt$p.value))

## clusterProfiler: GO and KEGG enrichment

cat("\nRunning clusterProfiler enrichment...\n")

# Helper: symbol → Entrez ID
symbols_to_entrez <- function(genes) {
  ids <- mapIds(org.Hs.eg.db, keys = genes,
                column = "ENTREZID", keytype = "SYMBOL",
                multiVals = "first")
  ids[!is.na(ids)]
}

run_enrichment <- function(gene_symbols, label, background = NULL) {
  entrez <- symbols_to_entrez(gene_symbols)
  if (length(entrez) < 5) {
    message("  ", label, ": too few genes (", length(entrez), ") — skipping")
    return(NULL)
  }

  bg_entrez <- if (!is.null(background)) symbols_to_entrez(background) else NULL

  # GO Biological Process
  go_res <- tryCatch(
    enrichGO(gene          = entrez,
             universe      = bg_entrez,
             OrgDb         = org.Hs.eg.db,
             ont           = "BP",
             pAdjustMethod = "fdr",
             pvalueCutoff  = 0.05,
             readable      = TRUE),
    error = function(e) { message("  GO failed: ", e$message); NULL }
  )

  # KEGG
  kegg_res <- tryCatch(
    enrichKEGG(gene          = entrez,
               universe      = bg_entrez,
               organism      = "hsa",
               pAdjustMethod = "fdr",
               pvalueCutoff  = 0.05),
    error = function(e) { message("  KEGG failed: ", e$message); NULL }
  )

  list(go = go_res, kegg = kegg_res, label = label,
       n_genes = length(entrez))
}

# Background: all intergenic VMR target genes (pooled)
all_intergenic_genes <- abc_links |>
  filter(tissue == "Pooled") |>
  pull(TargetGene_name) |>
  unique()

comparisons <- list(
  Pooled_Heritable    = gene_lists$Pooled$Heritable,
  Pooled_NonHeritable = gene_lists$Pooled$NonHeritable,
  Pooled_Q5           = genes_q5,
  Pooled_Q1Q4         = genes_q1q4
)

enrich_results <- lapply(names(comparisons), function(nm) {
  cat("  Running:", nm, "\n")
  run_enrichment(comparisons[[nm]], nm, background = all_intergenic_genes)
})
names(enrich_results) <- names(comparisons)

# Save enrichment tables
for (nm in names(enrich_results)) {
  er <- enrich_results[[nm]]
  if (is.null(er)) next

  nm_safe <- gsub("[^A-Za-z0-9]", "_", nm)

  if (!is.null(er$go)) {
    go_df <- as.data.frame(er$go)
    if (nrow(go_df) > 0) {
      go_df$comparison <- nm
      fwrite(go_df, file.path(OUT_DIR, paste0("enrichr_go_bp_", nm_safe, ".tsv")),
             sep = "\t")
    }
  }

  if (!is.null(er$kegg)) {
    kegg_df <- as.data.frame(er$kegg)
    if (nrow(kegg_df) > 0) {
      kegg_df$comparison <- nm
      fwrite(kegg_df, file.path(OUT_DIR, paste0("enrichr_kegg_", nm_safe, ".tsv")),
             sep = "\t")
    }
  }
}

# Combined GO results table (all comparisons)
go_all <- lapply(names(enrich_results), function(nm) {
  er <- enrich_results[[nm]]
  if (is.null(er) || is.null(er$go)) return(NULL)
  df <- as.data.frame(er$go)
  if (nrow(df) == 0) return(NULL)
  df$comparison <- nm
  df
}) |> bind_rows()

if (nrow(go_all) > 0) {
  fwrite(go_all, file.path(OUT_DIR, "enrichr_go_bp_all.tsv"), sep = "\t")
  cat("Saved: enrichr_go_bp_all.tsv\n")
}

kegg_all <- lapply(names(enrich_results), function(nm) {
  er <- enrich_results[[nm]]
  if (is.null(er) || is.null(er$kegg)) return(NULL)
  df <- as.data.frame(er$kegg)
  if (nrow(df) == 0) return(NULL)
  df$comparison <- nm
  df
}) |> bind_rows()

if (nrow(kegg_all) > 0) {
  fwrite(kegg_all, file.path(OUT_DIR, "enrichr_kegg_all.tsv"), sep = "\t")
  cat("Saved: enrichr_kegg_all.tsv\n")
}

cat("\nDone.\n")

#### Reproducibility ####
cat("\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
