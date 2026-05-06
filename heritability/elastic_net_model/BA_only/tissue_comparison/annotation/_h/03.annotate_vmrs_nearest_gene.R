#!/usr/bin/env Rscript

#### Nearest gene within 250 kb of VMRs (hg38 knownGene) ####

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(GenomicFeatures)
  library(here)

  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(sort(x), collapse = ";")
}

format_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x), x, paste0("chr", x))
}

make_region_key_gr <- function(gr) {
  paste(
    as.character(seqnames(gr)),
    start(gr),
    end(gr),
    sep = ":"
  )
}

load_vmrs <- function(tissue) {
  vmr_files <- here(
    paste0(
      "vmr-analysis/", tissue,
      "/_m/vmr/chr_",
      c(1:22, "X", "Y"),
      "/vmr.bed"
    )
  )

  vmr_list <- lapply(vmr_files, function(f) {
    if (file.exists(f)) {
      fread(f, col.names = c("chr", "start", "end"))
    } else {
      message("File does not exist: ", f)
      return(NULL)
    }
  })

  vmr_list <- Filter(Negate(is.null), vmr_list)

  if (length(vmr_list) == 0) {
    stop("No VMR files found for tissue: ", tissue)
  }

  vmr <- rbindlist(vmr_list)
  vmr[, chr := format_chr(chr)]

  vmr_gr <- GRanges(
    seqnames = vmr$chr,
    ranges = IRanges(start = vmr$start, end = vmr$end)
  )

  mcols(vmr_gr)$vmr_id <- paste0(tissue, "_vmr_", seq_along(vmr_gr))
  mcols(vmr_gr)$region_key <- make_region_key_gr(vmr_gr)

  return(vmr_gr)
}

load_gene_ranges <- function() {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  gene_gr <- GenomicFeatures::genes(txdb)
  gene_ids <- names(gene_gr)

  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_ids,
    keytype = "ENTREZID",
    columns = c("SYMBOL", "GENENAME")
  ) %>%
    distinct(ENTREZID, .keep_all = TRUE)

  mcols(gene_gr)$gene_id <- gene_ids
  mcols(gene_gr)$gene_symbol <- gene_map$SYMBOL[
    match(gene_ids, gene_map$ENTREZID)
  ]
  mcols(gene_gr)$gene_name <- gene_map$GENENAME[
    match(gene_ids, gene_map$ENTREZID)
  ]

  mcols(gene_gr)$gene_symbol <- ifelse(
    is.na(mcols(gene_gr)$gene_symbol),
    mcols(gene_gr)$gene_id,
    mcols(gene_gr)$gene_symbol
  )

  GenomeInfoDb::seqlevelsStyle(gene_gr) <- "UCSC"

  return(gene_gr)
}

prepare_enet <- function(enet) {
  enet <- na.omit(enet)

  enet <- enet %>%
    dplyr::select(chrom, start, end, h2_unscaled, r_squared_cv) %>%
    mutate(
      chr = format_chr(chrom),
      region_key = paste(chr, start, end, sep = ":"),
      h2_category = case_when(
        r_squared_cv < 0.3 ~ "Low prediction",
        h2_unscaled < 0.1 & r_squared_cv >= 0.3 ~ "Non-heritable",
        h2_unscaled >= 0.1 & r_squared_cv >= 0.3 ~ "Heritable"
      ),
      h2_category = factor(
        h2_category,
        levels = c("Heritable", "Non-heritable", "Low prediction")
      )
    ) %>%
    dplyr::select(
      region_key,
      h2_unscaled,
      r_squared_cv,
      h2_category
    )

  return(enet)
}

annotate_genes_within_window <- function(
  vmr_gr,
  gene_gr,
  window_bp = 250000
) {
  common_seqlevels <- intersect(seqlevels(vmr_gr), seqlevels(gene_gr))

  vmr_gr <- keepSeqlevels(
    vmr_gr,
    common_seqlevels,
    pruning.mode = "coarse"
  )

  gene_gr <- keepSeqlevels(
    gene_gr,
    common_seqlevels,
    pruning.mode = "coarse"
  )

  vmr_window_gr <- GRanges(
    seqnames = seqnames(vmr_gr),
    ranges = IRanges(
      start = pmax(1, start(vmr_gr) - window_bp),
      end = end(vmr_gr) + window_bp
    ),
    strand = strand(vmr_gr)
  )

  mcols(vmr_window_gr)$vmr_id <- mcols(vmr_gr)$vmr_id
  mcols(vmr_window_gr)$region_key <- mcols(vmr_gr)$region_key

  hits <- findOverlaps(
    vmr_window_gr,
    gene_gr,
    ignore.strand = TRUE
  )

  base_df <- data.frame(
    seqnames = as.character(seqnames(vmr_gr)),
    start = start(vmr_gr),
    end = end(vmr_gr),
    width = width(vmr_gr),
    strand = as.character(strand(vmr_gr)),
    vmr_id = mcols(vmr_gr)$vmr_id,
    region_key = mcols(vmr_gr)$region_key
  )

  if (length(hits) == 0) {
    out <- base_df %>%
      mutate(
        n_genes_within_250kb = 0L,
        gene_ids_within_250kb = NA_character_,
        gene_symbols_within_250kb = NA_character_,
        gene_names_within_250kb = NA_character_,
        gene_distances_within_250kb = NA_character_,
        nearest_gene_id_within_250kb = NA_character_,
        nearest_gene_symbol_within_250kb = NA_character_,
        nearest_gene_name_within_250kb = NA_character_,
        distance_to_nearest_gene_within_250kb = NA_integer_
      )

    return(out)
  }

  qh <- queryHits(hits)
  sh <- subjectHits(hits)

  pair_df <- data.table(
    region_key = mcols(vmr_gr)$region_key[qh],
    gene_id = mcols(gene_gr)$gene_id[sh],
    gene_symbol = mcols(gene_gr)$gene_symbol[sh],
    gene_name = mcols(gene_gr)$gene_name[sh],
    gene_chr = as.character(seqnames(gene_gr))[sh],
    gene_start = start(gene_gr)[sh],
    gene_end = end(gene_gr)[sh],
    gene_strand = as.character(strand(gene_gr))[sh],
    distance_to_gene = GenomicRanges::distance(
      vmr_gr[qh], gene_gr[sh], ignore.strand = TRUE
    )
  )

  pair_df <- pair_df[distance_to_gene <= window_bp]
  pair_df <- pair_df[order(region_key, distance_to_gene, gene_symbol)]

  nearest_df <- pair_df[
    ,
    .SD[1],
    by = region_key
  ][
    ,
    .(
      region_key,
      nearest_gene_id_within_250kb = gene_id,
      nearest_gene_symbol_within_250kb = gene_symbol,
      nearest_gene_name_within_250kb = gene_name,
      nearest_gene_chr_within_250kb = gene_chr,
      nearest_gene_start_within_250kb = gene_start,
      nearest_gene_end_within_250kb = gene_end,
      nearest_gene_strand_within_250kb = gene_strand,
      distance_to_nearest_gene_within_250kb = distance_to_gene
    )
  ]

  collapsed_df <- pair_df[
    ,
    .(
      n_genes_within_250kb = uniqueN(gene_id),
      gene_ids_within_250kb = collapse_unique(gene_id),
      gene_symbols_within_250kb = collapse_unique(gene_symbol),
      gene_names_within_250kb = collapse_unique(gene_name),
      gene_distances_within_250kb = paste(
        paste0(gene_symbol, ":", distance_to_gene),
        collapse = ";"
      )
    ),
    by = region_key
  ]

  out <- base_df %>%
    left_join(collapsed_df, by = "region_key") %>%
    left_join(nearest_df, by = "region_key") %>%
    mutate(
      n_genes_within_250kb = ifelse(
        is.na(n_genes_within_250kb),
        0L,
        n_genes_within_250kb
      )
    )

  return(out)
}

annotate_vmrs_with_genes_250kb <- function(
  vmr_gr, enet, gene_gr, out_file, window_bp = 250000
) {
  gene_annot <- annotate_genes_within_window(
    vmr_gr = vmr_gr, gene_gr = gene_gr, window_bp = window_bp
  )

  enet <- prepare_enet(enet) |>
    left_join(gene_annot, by = "region_key") |>
    arrange(seqnames, start, end)

  fwrite(enet, out_file, sep = "\t")

  return(enet)
}

out_path <- here(
  "heritability/elastic_net_model/BA_only/tissue_comparison/annotation/_m"
)

if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

tissues <- c("dlpfc", "caudate", "hippocampus")
window_bp <- 250000
gene_gr <- load_gene_ranges()

for (tissue in tissues) {
  message("Processing tissue: ", tissue)

  out_gene <- file.path(
    out_path,
    paste0(tissue, "_vmr_genes_within_250kb_hg38.tsv")
  )

  vmr_gr <- load_vmrs(tissue)

  enet_file <- here(
    "heritability/elastic_net_model/BA_only",
    paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv")
  )

  enet <- fread(enet_file)

  vmr_gene_annot <- annotate_vmrs_with_genes_250kb(
    vmr_gr = vmr_gr,
    enet = enet,
    gene_gr = gene_gr,
    out_file = out_gene,
    window_bp = window_bp
  )

  summary_df <- vmr_gene_annot %>%
    group_by(h2_category) %>%
    summarise(
      n_vmrs = n(),
      n_vmrs_with_gene_within_250kb = sum(n_genes_within_250kb > 0),
      n_vmrs_without_gene_within_250kb = sum(n_genes_within_250kb == 0),
      median_n_genes_within_250kb = median(n_genes_within_250kb),
      .groups = "drop"
    )

  print(summary_df)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
if (requireNamespace("sessioninfo", quietly = TRUE)) {
  sessioninfo::session_info()
}
