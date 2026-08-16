#### Intergenic VMR ~ Repressive Chromatin Overlap Analysis ####
##
## Tests whether heritable VMRs — especially heritable intergenic VMRs — are
## enriched in repressive chromatin domains relative to non-heritable VMRs.
##
## Repressive features tested:
##   H3K27me3  (Polycomb silencing)       — Roadmap gappedPeak (primary)
##   H3K9me3   (constitutive heterochrom) — Roadmap gappedPeak (primary)
##   ChromHMM 15-state repressive states:
##     13_ReprPC + 14_ReprPCWk  → Polycomb-repressed
##     9_Het                    → Heterochromatin
##     15_Quies                 → Quiescent/Low
##     union of above           → BroadRepressive
##     10_TssBiv + 11_BivFlnk + 12_EnhBiv → Bivalent (reported separately)
##
## Comparisons:
##   1. All heritable vs all non-heritable VMRs (Fisher's exact)
##   2. Heritable intergenic vs non-heritable intergenic VMRs (Fisher's exact)
##   3. Continuous h² ~ repressive overlap (logistic regression + covariates)
##   4. Top h² quintile (Q5) vs Q1-Q4 within intergenic VMRs
##   5. Top h² decile vs bottom 90% within intergenic VMRs
##
## Sensitivity:
##   broadPeak vs gappedPeak for H3K27me3 and H3K9me3
##
## Tissue–Epigenome mapping (NIH Roadmap consolidated epigenomes):
##   Caudate     → E068 (Anterior Caudate)
##   DLPFC       → E073 (Dorsolateral Prefrontal Cortex)
##   Hippocampus → E071 (Hippocampus Middle)
##
## Coordinate liftover: VMRs in hg38 → hg19 using rtracklayer::liftOver
## Chain: inputs/supportfiles/_m/hg38ToHg19.over.chain

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

## Configuration

TISSUES     <- c("Caudate", "DLPFC", "Hippocampus")
R2_THRESH   <- 0.3
N_QUINTILES <- 5
H2_UNIT     <- 0.1

# Tissue → Roadmap EID mapping
EID_MAP <- c(
  Caudate     = "E068",   # Anterior Caudate
  DLPFC       = "E073",   # Dorsolateral Prefrontal Cortex
  Hippocampus = "E071"    # Hippocampus Middle
)

# ChromHMM state groupings (15-state core model)
CHROMHMM_GROUPS <- list(
  Polycomb     = c("13_ReprPC", "14_ReprPCWk"),
  Het          = c("9_Het"),
  Quies        = c("15_Quies"),
  Bivalent     = c("10_TssBiv", "11_BivFlnk", "12_EnhBiv"),
  BroadRepressive = c("9_Het", "13_ReprPC", "14_ReprPCWk", "15_Quies")
)

## Paths

ANNOT_DIR   <- here("heritability", "elastic_net_model", "BA_only",
                    "tissue_comparison", "annotation", "_m")
ENET_BASE   <- here("heritability", "elastic_net_model", "BA_only")
CHAIN_FILE  <- here("inputs", "supportfiles", "_m", "hg38ToHg19.over.chain")
ROADMAP_DIR <- file.path( ## Move later to inputs/supportfiles/_m
  "/projects/b1213/resources/EpigenomeRoadmap",
  "egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated"
)
CHROMHMM_DIR <- "/projects/b1213/resources/EpigenomeRoadmap/hmm/15core"
OUT_DIR      <- here("heritability", "elastic_net_model", "BA_only",
                     "tissue_comparison", "annotation", "repressive_chromatin", "_m")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

chain <- import.chain(CHAIN_FILE)

## Functions

safe_glm <- function(formula, data) {
  tryCatch(
    glm(formula, data = data, family = binomial),
    error = function(e) { warning("GLM failed: ", conditionMessage(e)); NULL }
  )
}

# Agresti-Coull 95% CI
ac_ci <- function(n_annot, n, z = qnorm(0.975)) {
  n_tilde <- n + z^2
  p_tilde <- (n_annot + z^2 / 2) / n_tilde
  lo <- pmax(0, p_tilde - z * sqrt(p_tilde * (1 - p_tilde) / n_tilde))
  hi <- pmin(1, p_tilde + z * sqrt(p_tilde * (1 - p_tilde) / n_tilde))
  list(lo = lo, hi = hi)
}

## Data loading

load_vmrs <- function(tissue, intergenic_only = FALSE) {
  tissue_lower <- tolower(tissue)

  annot <- fread(file.path(ANNOT_DIR,
    paste0(tissue_lower, "_vmr_annotations_hg38_wide.tsv")))

  enet <- fread(file.path(ENET_BASE, tissue_lower, "_m",
    paste0(tissue_lower, "_summary_elastic-net.tsv"))) |>
    mutate(seqnames = paste0("chr", chrom)) |>
    dplyr::select(seqnames, start, end, h2_unscaled, r_squared_cv, num_snps)

  df <- annot |>
    left_join(enet, by = c("seqnames", "start", "end")) |>
    mutate(
      vmr_length = end - start,
      h2_scaled  = h2_unscaled / H2_UNIT,
      tissue     = tissue
    ) |>
    filter(
      !is.na(r_squared_cv), r_squared_cv > R2_THRESH,
      !is.na(h2_unscaled), is.finite(h2_unscaled),
      vmr_length > 0, !is.na(num_snps), num_snps > 0,
      h2_category %in% c("Heritable", "Non-heritable")
    )

  if (intergenic_only) {
    df <- filter(df, hg38_genes_intergenic == 1)
  }
  df
}

liftover_df <- function(df) {
  gr38 <- GRanges(
    seqnames = df$seqnames,
    ranges   = IRanges(start = df$start, end = df$end),
    idx      = seq_len(nrow(df))
  )
  lo <- liftOver(gr38, chain)
  # Keep only VMRs that lifted to a single range (exclude multi-mappers / failures)
  n_ranges <- lengths(lo)
  keep_idx <- which(n_ranges == 1)
  gr19 <- unlist(lo[keep_idx])
  n_dropped <- nrow(df) - length(keep_idx)
  if (n_dropped > 0) {
    cat(sprintf("    liftOver: dropped %d VMRs (%.1f%%) that did not map uniquely to hg19\n",
                n_dropped, 100 * n_dropped / nrow(df)))
  }
  df_kept <- df[keep_idx, ]
  df_kept$seqnames_hg19 <- as.character(seqnames(gr19))
  df_kept$start_hg19    <- start(gr19)
  df_kept$end_hg19      <- end(gr19)
  df_kept
}

load_roadmap_peaks <- function(eid, mark, peak_type = "gappedPeak") {
  fname <- sprintf("%s-%s.%s.gz", eid, mark, peak_type)
  fpath <- file.path(ROADMAP_DIR, peak_type, fname)
  if (!file.exists(fpath)) {
    stop("Peak file not found: ", fpath)
  }
  # gappedPeak is BED12-like; broadPeak is BED6+4
  # Both have chr/start/end in cols 1-3 — read only those
  dt <- fread(cmd = paste("zcat", shQuote(fpath)), header = FALSE,
              select = 1:3, col.names = c("chr", "start", "end"))
  GRanges(seqnames = dt$chr, ranges = IRanges(start = dt$start, end = dt$end))
}

load_chromhmm <- function(eid, states) {
  fpath <- file.path(CHROMHMM_DIR, sprintf("%s_15_coreMarks_mnemonics.bed.gz", eid))
  if (!file.exists(fpath)) stop("ChromHMM file not found: ", fpath)
  dt <- fread(cmd = paste("zcat", shQuote(fpath)), header = FALSE,
              col.names = c("chr", "start", "end", "state"))
  dt_filt <- dt[state %in% states]
  GRanges(seqnames = dt_filt$chr,
          ranges   = IRanges(start = dt_filt$start, end = dt_filt$end))
}

flag_overlaps <- function(df, annot_gr_list) {
  # Build hg19 GRanges from df
  vmr_gr <- GRanges(
    seqnames = df$seqnames_hg19,
    ranges   = IRanges(start = df$start_hg19, end = df$end_hg19)
  )
  for (ann_name in names(annot_gr_list)) {
    hits    <- findOverlaps(vmr_gr, annot_gr_list[[ann_name]], ignore.strand = TRUE)
    in_peak <- logical(nrow(df))
    in_peak[queryHits(hits)] <- TRUE
    df[[paste0("in_", ann_name)]] <- as.integer(in_peak)
  }
  df
}

run_fishers <- function(df, ann_col, tissue, comparison) {
  vals   <- df[[ann_col]]
  h2_cat <- df$h2_category
  tab    <- table(vals, h2_cat)
  # Ensure both levels present
  if (!all(c("Heritable", "Non-heritable") %in% colnames(tab))) return(NULL)
  if (nrow(tab) < 2) return(NULL)
  tab <- tab[, c("Heritable", "Non-heritable"), drop = FALSE]
  ft  <- fisher.test(tab)
  tibble(
    tissue      = tissue,
    comparison  = comparison,
    annotation  = gsub("^in_", "", ann_col),
    or          = ft$estimate,
    ci_lo       = ft$conf.int[1],
    ci_hi       = ft$conf.int[2],
    p_value     = ft$p.value,
    n_heritable = sum(h2_cat == "Heritable"),
    n_nonher    = sum(h2_cat == "Non-heritable"),
    n_in_her    = sum(vals == 1 & h2_cat == "Heritable"),
    n_in_nonher = sum(vals == 1 & h2_cat == "Non-heritable")
  )
}

run_logistic <- function(df, ann_col, tissue, comparison) {
  df_mod <- df |>
    mutate(
      y        = .data[[ann_col]],
      log_len  = log(vmr_length),
      log_snps = log(num_snps + 1)
    )
  fit <- safe_glm(y ~ h2_unscaled + log_len + log_snps, data = df_mod)
  if (is.null(fit)) return(NULL)
  coef_df <- as.data.frame(summary(fit)$coefficients)
  coef_df$term <- rownames(coef_df); rownames(coef_df) <- NULL
  ci <- tryCatch(confint(fit), error = function(e) NULL)
  h2_row <- coef_df[coef_df$term == "h2_unscaled", ]
  if (nrow(h2_row) == 0) return(NULL)
  tibble(
    tissue      = tissue,
    comparison  = comparison,
    annotation  = gsub("^in_", "", ann_col),
    estimate    = h2_row$Estimate,
    std_error   = h2_row$`Std. Error`,
    z_value     = h2_row$`z value`,
    p_value     = h2_row$`Pr(>|z|)`,
    or          = exp(h2_row$Estimate),
    ci_lo_log   = if (!is.null(ci)) ci["h2_unscaled", 1] else NA_real_,
    ci_hi_log   = if (!is.null(ci)) ci["h2_unscaled", 2] else NA_real_,
    ci_lo_or    = if (!is.null(ci)) exp(ci["h2_unscaled", 1]) else NA_real_,
    ci_hi_or    = if (!is.null(ci)) exp(ci["h2_unscaled", 2]) else NA_real_,
    n_vmrs      = nrow(df_mod),
    n_in_peak   = sum(df_mod$y, na.rm = TRUE)
  )
}

run_topbin_fishers <- function(df, ann_col, tissue, comparison, top_frac) {
  thresh <- quantile(df$h2_unscaled, probs = 1 - top_frac, na.rm = TRUE)
  df2    <- df |>
    mutate(is_top = as.integer(h2_unscaled > thresh))
  tab <- table(df2[[ann_col]], df2$is_top)
  if (!all(c("0", "1") %in% colnames(tab))) return(NULL)
  if (nrow(tab) < 2) return(NULL)
  tab <- tab[, c("0", "1"), drop = FALSE]
  ft  <- fisher.test(tab)
  label <- if (top_frac == 0.2) "Q5_vs_Q1Q4" else "Decile_vs_Rest"
  tibble(
    tissue      = tissue,
    comparison  = paste0(comparison, "_", label),
    annotation  = gsub("^in_", "", ann_col),
    or          = ft$estimate,
    ci_lo       = ft$conf.int[1],
    ci_hi       = ft$conf.int[2],
    p_value     = ft$p.value,
    n_top       = sum(df2$is_top == 1),
    n_lower     = sum(df2$is_top == 0),
    n_in_top    = sum(df2[[ann_col]] == 1 & df2$is_top == 1),
    n_in_lower  = sum(df2[[ann_col]] == 1 & df2$is_top == 0)
  )
}

chromhmm_rate_table <- function(df, tissue, all_states, eid) {
  # Load all 15 states
  fpath <- file.path(CHROMHMM_DIR, sprintf("%s_15_coreMarks_mnemonics.bed.gz", eid))
  dt    <- fread(cmd = paste("zcat", shQuote(fpath)), header = FALSE,
                 col.names = c("chr", "start", "end", "state"))
  vmr_gr <- GRanges(seqnames = df$seqnames_hg19,
                    ranges   = IRanges(start = df$start_hg19, end = df$end_hg19))
  rows <- lapply(unique(dt$state), function(st) {
    gr_st   <- GRanges(seqnames = dt[state == st, chr],
                       ranges   = IRanges(dt[state == st, start], dt[state == st, end]))
    hits    <- findOverlaps(vmr_gr, gr_st, ignore.strand = TRUE)
    in_st   <- logical(nrow(df))
    in_st[queryHits(hits)] <- TRUE
    df$in_state <- as.integer(in_st)
    df |>
      group_by(h2_category) |>
      summarise(
        n        = n(),
        n_in     = sum(in_state),
        frac     = n_in / n,
        .groups  = "drop"
      ) |>
      mutate(state = st, tissue = tissue)
  })
  bind_rows(rows)
}

## Main analysis

cat("Loading Roadmap peak GRanges...\n")
peaks_gr <- list()
for (tissue in TISSUES) {
  eid <- EID_MAP[tissue]
  cat(sprintf("  %s (%s)\n", tissue, eid))
  for (mark in c("H3K27me3", "H3K9me3")) {
    for (pt in c("gappedPeak", "broadPeak")) {
      key <- paste(tissue, mark, pt, sep = "_")
      peaks_gr[[key]] <- load_roadmap_peaks(eid, mark, pt)
      cat(sprintf("    %s %s (%s): %d peaks\n",
                  mark, pt, eid, length(peaks_gr[[key]])))
    }
  }
}

cat("\nLoading ChromHMM GRanges...\n")
chromhmm_gr <- list()
for (tissue in TISSUES) {
  eid <- EID_MAP[tissue]
  for (grp in names(CHROMHMM_GROUPS)) {
    key <- paste(tissue, grp, sep = "_")
    chromhmm_gr[[key]] <- load_chromhmm(eid, CHROMHMM_GROUPS[[grp]])
    cat(sprintf("  %s %s (%s): %d regions\n",
                tissue, grp, eid, length(chromhmm_gr[[key]])))
  }
}

cat("\nLoading VMR data and lifting over coordinates...\n")
vmr_all_list   <- list()
vmr_inter_list <- list()
for (tissue in TISSUES) {
  cat(sprintf("  %s\n", tissue))
  df_all   <- load_vmrs(tissue, intergenic_only = FALSE)
  df_inter <- df_all |> filter(hg38_genes_intergenic == 1)

  cat(sprintf("    All: %d VMRs (%d Her, %d Non-her) before liftover\n",
              nrow(df_all), sum(df_all$h2_category == "Heritable"),
              sum(df_all$h2_category == "Non-heritable")))

  df_all   <- liftover_df(df_all)
  df_inter <- liftover_df(df_inter)

  cat(sprintf("    Intergenic after liftover: %d (%d Her, %d Non-her)\n",
              nrow(df_inter), sum(df_inter$h2_category == "Heritable"),
              sum(df_inter$h2_category == "Non-heritable")))

  # Build per-tissue annotation GRanges list
  eid <- EID_MAP[tissue]
  ann_list <- list(
    H3K27me3           = peaks_gr[[paste(tissue, "H3K27me3", "gappedPeak", sep = "_")]],
    H3K9me3            = peaks_gr[[paste(tissue, "H3K9me3",  "gappedPeak", sep = "_")]],
    H3K27me3_broad     = peaks_gr[[paste(tissue, "H3K27me3", "broadPeak",  sep = "_")]],
    H3K9me3_broad      = peaks_gr[[paste(tissue, "H3K9me3",  "broadPeak",  sep = "_")]],
    Polycomb           = chromhmm_gr[[paste(tissue, "Polycomb",        sep = "_")]],
    Het                = chromhmm_gr[[paste(tissue, "Het",             sep = "_")]],
    Quies              = chromhmm_gr[[paste(tissue, "Quies",           sep = "_")]],
    BroadRepressive    = chromhmm_gr[[paste(tissue, "BroadRepressive", sep = "_")]],
    Bivalent           = chromhmm_gr[[paste(tissue, "Bivalent",        sep = "_")]]
  )

  df_all   <- flag_overlaps(df_all,   ann_list)
  df_inter <- flag_overlaps(df_inter, ann_list)

  vmr_all_list[[tissue]]   <- df_all
  vmr_inter_list[[tissue]] <- df_inter
}

## Save annotated VMR table
cat("\nSaving annotated VMR table...\n")
vmr_all_out <- bind_rows(vmr_all_list)
fwrite(vmr_all_out, file.path(OUT_DIR, "vmr_repressive_overlap.tsv"), sep = "\t")
cat("Saved: vmr_repressive_overlap.tsv\n")

# Run enrichment tests

# Annotation columns to test
ANN_COLS_PEAK  <- c("in_H3K27me3", "in_H3K9me3",
                    "in_H3K27me3_broad", "in_H3K9me3_broad")
ANN_COLS_HMM   <- c("in_Polycomb", "in_Het", "in_Quies",
                    "in_BroadRepressive", "in_Bivalent")
ANN_COLS_ALL   <- c(ANN_COLS_PEAK, ANN_COLS_HMM)

# Fisher's exact
cat("\nRunning Fisher's exact tests...\n")
fish_list <- list()

for (tissue in TISSUES) {
  # Comparison 1: all VMRs
  for (col in ANN_COLS_ALL) {
    res <- run_fishers(vmr_all_list[[tissue]], col, tissue, "All_VMRs")
    if (!is.null(res)) fish_list[[length(fish_list) + 1]] <- res
  }
  # Comparison 2: intergenic VMRs
  for (col in ANN_COLS_ALL) {
    res <- run_fishers(vmr_inter_list[[tissue]], col, tissue, "Intergenic_VMRs")
    if (!is.null(res)) fish_list[[length(fish_list) + 1]] <- res
  }
}

fisher_df <- bind_rows(fish_list)
fisher_df$fdr <- p.adjust(fisher_df$p_value, method = "fdr")
fwrite(fisher_df, file.path(OUT_DIR, "fishers_repressive_enrichment.tsv"), sep = "\t")
cat("Saved: fishers_repressive_enrichment.tsv\n")

# Logistic regression
cat("\nRunning logistic regressions (continuous h²)...\n")
logistic_list <- list()

for (tissue in TISSUES) {
  for (col in ANN_COLS_ALL) {
    res <- run_logistic(vmr_all_list[[tissue]],   col, tissue, "All_VMRs")
    if (!is.null(res)) logistic_list[[length(logistic_list) + 1]] <- res
    res <- run_logistic(vmr_inter_list[[tissue]], col, tissue, "Intergenic_VMRs")
    if (!is.null(res)) logistic_list[[length(logistic_list) + 1]] <- res
  }
}

logistic_df <- bind_rows(logistic_list)
logistic_df$fdr <- p.adjust(logistic_df$p_value, method = "fdr")
fwrite(logistic_df, file.path(OUT_DIR, "logistic_repressive_enrichment.tsv"), sep = "\t")
cat("Saved: logistic_repressive_enrichment.tsv\n")

# Q5 and decile Fisher's exact
cat("\nRunning top-quintile and top-decile Fisher's tests (intergenic VMRs)...\n")
topbin_list <- list()

for (tissue in TISSUES) {
  df_inter <- vmr_inter_list[[tissue]]
  for (col in ANN_COLS_ALL) {
    # Q5 vs Q1-Q4
    res <- run_topbin_fishers(df_inter, col, tissue, "Intergenic_VMRs", top_frac = 0.2)
    if (!is.null(res)) topbin_list[[length(topbin_list) + 1]] <- res
    # Top decile vs rest
    res <- run_topbin_fishers(df_inter, col, tissue, "Intergenic_VMRs", top_frac = 0.1)
    if (!is.null(res)) topbin_list[[length(topbin_list) + 1]] <- res
  }
}

topbin_df <- bind_rows(topbin_list)
topbin_df$fdr <- p.adjust(topbin_df$p_value, method = "fdr")
fwrite(topbin_df, file.path(OUT_DIR, "q5_decile_repressive_fishers.tsv"), sep = "\t")
cat("Saved: q5_decile_repressive_fishers.tsv\n")

# ChromHMM state overlap rate table
cat("\nBuilding ChromHMM state overlap rate table (all VMRs)...\n")
rate_list <- lapply(TISSUES, function(tissue) {
  eid <- EID_MAP[tissue]
  chromhmm_rate_table(vmr_all_list[[tissue]], tissue, CHROMHMM_GROUPS, eid)
})
rate_df <- bind_rows(rate_list)
fwrite(rate_df, file.path(OUT_DIR, "chromhmm_state_overlap_rates.tsv"), sep = "\t")
cat("Saved: chromhmm_state_overlap_rates.tsv\n")

# Summary print

cat("\n=== Fisher's enrichment: Intergenic heritable vs non-heritable ===\n")
fisher_df |>
  filter(comparison == "Intergenic_VMRs",
         !grepl("_broad", annotation),
         annotation != "H3K27me3_broad", annotation != "H3K9me3_broad") |>
  dplyr::select(tissue, annotation, or, ci_lo, ci_hi, p_value, fdr,
                n_in_her, n_in_nonher) |>
  arrange(fdr) |>
  print(n = Inf)

cat("\n=== Logistic regression: continuous h² ~ repressive overlap (intergenic) ===\n")
logistic_df |>
  filter(comparison == "Intergenic_VMRs") |>
  dplyr::select(tissue, annotation, or, ci_lo_or, ci_hi_or, p_value, fdr) |>
  arrange(fdr) |>
  print(n = Inf)

#### Reproducibility ####
cat("\nReproducibility information:\n")
print(Sys.time())
print(proc.time())
options(width = 120)
sessioninfo::session_info()
