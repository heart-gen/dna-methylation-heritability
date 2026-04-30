#### Intergenic VMR ~ Brain Open Chromatin Overlap Analysis ####
##
## Tests whether high-heritability intergenic VMRs are preferentially
## enriched in brain cell-type open chromatin regulatory elements (CREs).
##
## Comparisons:
##   1. Heritable vs Non-heritable intergenic VMRs (categorical; Fisher's exact)
##   2. Continuous h² ~ open chromatin overlap (logistic regression)
##   3. Top h² quintile (Q5) vs lower quintiles (Q1–Q4) within intergenic VMRs
##
## Open chromatin resources: BrainScope cell-type ATAC-seq peaks + union
## Cell types: Exc, Inh, Astro, Oligo, OPC, Endo, Micro (+ union)
##
## Run: conda run -p $ENV_PATH/epigenomics Rscript ../_h/01.overlap_open_chromatin.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(GenomicRanges)
})

## Configuration

TISSUES     <- c("Caudate", "DLPFC", "Hippocampus")
POPULATIONS <- c("AA", "EA")
R2_THRESH   <- 0.3
N_QUINTILES <- 5
H2_UNIT     <- 0.1

ANNOT_DIR <- here::here(
  "heritability", "elastic_net_model", "all_individuals",
  "tissue_comparison", "annotation", "_m"
)
ENET_BASE <- here::here(
  "heritability", "elastic_net_model", "all_individuals"
)
BRAINSCOPE_DIR <- here::here("inputs", "brainscope", "_m")
OUT_DIR <- here::here(
  "heritability", "elastic_net_model", "all_individuals",
  "tissue_comparison", "annotation", "open_chromatin", "_m"
)

# Cell-type ATAC peak files
CELL_TYPES <- c(
  "Union"  = "All.celltypes.Union.PeakCalls.bed",
  "Exc"    = "Exc.PeakCalls.bed",
  "Inh"    = "Inh.PeakCalls.bed",
  "Astro"  = "Astro.PeakCalls.bed",
  "Oligo"  = "Oligo.PeakCalls.bed",
  "OPC"    = "OPC.PeakCalls.bed",
  "Endo"   = "Endo.PeakCalls.bed",
  "Micro"  = "Micro.PeakCalls.bed"
)

## Helpers

load_atac_peaks <- function(bed_file) {
  df <- fread(file.path(BRAINSCOPE_DIR, bed_file), header = FALSE)
  # Handle variable column counts: 3 cols (chr, start, end) or
  # 4 cols (index, chr, start, end) — detect by whether col1 is numeric
  if (ncol(df) >= 4 && !grepl("^chr", df[[1]][1])) {
    # First column is an index
    chr_col   <- 2L
    start_col <- 3L
    end_col   <- 4L
  } else {
    chr_col   <- 1L
    start_col <- 2L
    end_col   <- 3L
  }
  GRanges(
    seqnames = df[[chr_col]],
    ranges   = IRanges(start = df[[start_col]], end = df[[end_col]])
  )
}

# Agresti-Coull 95% CI bounds
ac_ci_lo <- function(n_annot, n, z = qnorm(0.975)) {
  n_tilde <- n + z^2
  p_tilde <- (n_annot + z^2 / 2) / n_tilde
  pmax(0, p_tilde - z * sqrt(p_tilde * (1 - p_tilde) / n_tilde))
}
ac_ci_hi <- function(n_annot, n, z = qnorm(0.975)) {
  n_tilde <- n + z^2
  p_tilde <- (n_annot + z^2 / 2) / n_tilde
  pmin(1, p_tilde + z * sqrt(p_tilde * (1 - p_tilde) / n_tilde))
}

safe_glm <- function(formula, data) {
  tryCatch(
    glm(formula, data = data, family = binomial),
    error = function(e) {
      warning("GLM failed: ", conditionMessage(e))
      NULL
    }
  )
}

load_intergenic_vmrs <- function(tissue, pop) {
  tissue_lower <- tolower(tissue)

  # Wide annotations: binary annotation columns, one row per VMR
  annot <- fread(file.path(ANNOT_DIR,
    paste0(tissue_lower, "_vmr_annotations_hg38_wide.tsv")))

  # Elastic-net: h2_unscaled, r_squared_cv, num_snps
  enet <- fread(file.path(ENET_BASE, tissue_lower, "_m",
    paste0(tissue_lower, "_summary_elastic-net_matched_r2_0.3.tsv"))) |>
    dplyr::rename(seqnames = chrom) |>
    dplyr::select(seqnames, start, end,
                  h2_unscaled = paste0("h2_unscaled_", pop),
                  r_squared_cv = paste0("r_squared_cv_", pop),
                  num_snps = paste0("num_snps_", pop))

  df <- annot |>
    left_join(enet, by = c("seqnames", "start", "end")) |>
    mutate(
      vmr_length = end - start,
      h2_scaled  = h2_unscaled / H2_UNIT,
      tissue     = tissue,
      population = pop
    )

  # Filter to intergenic VMRs with high-confidence h2 estimates
  df |>
    filter(
      hg38_genes_intergenic == 1,
      !is.na(r_squared_cv), r_squared_cv > R2_THRESH,
      !is.na(h2_unscaled), is.finite(h2_unscaled),
      vmr_length > 0, !is.na(num_snps), num_snps > 0,
      h2_category %in% c("Heritable", "Non-heritable")
    )
}


## Load ATAC peaks once

cat("\nLoading cell-type ATAC peaks...\n")
atac_gr <- lapply(CELL_TYPES, function(fn) {
  gr <- load_atac_peaks(fn)
  cat(sprintf("  %-8s: %d peaks\n", names(CELL_TYPES)[CELL_TYPES == fn], length(gr)))
  gr
})

names(atac_gr) <- names(CELL_TYPES)

for (pop in POPULATIONS) {

  cat(sprintf("Processing population: %s \n", pop))

  ## Load all tissues

  vmr_list <- lapply(TISSUES, function(t) load_intergenic_vmrs(t, pop)) 

  names(vmr_list) <- TISSUES 

  for (t in TISSUES) { 
    df <- vmr_list[[t]] 
    cat(sprintf(" %s: %d intergenic VMRs (%d Heritable, %d Non-heritable)\n", t, 
    nrow(df), 
    sum(df$h2_category == "Heritable"),
     sum(df$h2_category == "Non-heritable") )) 
  }

  ## Overlap: flag each VMR as overlapping each cell type's peaks

  cat("\nComputing overlaps...\n")

  vmr_annotated_list <- lapply(TISSUES, function(t) {
    df <- vmr_list[[t]]
    vmr_gr <- GRanges(
      seqnames = df$seqnames,
      ranges   = IRanges(start = df$start, end = df$end)
    )

    for (ct in names(atac_gr)) {
      hits <- findOverlaps(vmr_gr, atac_gr[[ct]], ignore.strand = TRUE)
      in_peak <- logical(nrow(df))
      in_peak[queryHits(hits)] <- TRUE
      df[[paste0("in_", ct)]] <- as.integer(in_peak)
    }

    # Number of cell types (excluding Union) with overlap
    ct_cols <- paste0("in_", setdiff(names(atac_gr), "Union"))
    df$n_celltypes_overlap <- rowSums(df[, ct_cols, with = FALSE])

    df
  })
  names(vmr_annotated_list) <- TISSUES

  # Save annotated intergenic VMR table
  vmr_all <- bind_rows(vmr_annotated_list)
  fwrite(vmr_all, file.path(OUT_DIR, paste0("intergenic_vmr_atac_overlap_", pop, ".tsv")), sep = "\t")
  cat(sprintf("Saved: intergenic_vmr_atac_overlap_%s.tsv\n", pop))

  ## Fisher's exact test: Heritable vs Non-heritable per cell type × tissue

  cat("\nRunning Fisher's exact tests (per tissue)...\n")

  run_fishers <- function(df, ct_col, tissue) {
    ct_vals   <- df[[ct_col]]
    h2_cat    <- df$h2_category

    tab <- table(ct_vals, h2_cat)[, c("Heritable", "Non-heritable")]
    if (nrow(tab) < 2) return(NULL)

    ft <- fisher.test(tab)
    tibble(
      tissue      = tissue,
      cell_type   = gsub("^in_", "", ct_col),
      or          = ft$estimate,
      ci_lo       = ft$conf.int[1],
      ci_hi       = ft$conf.int[2],
      p_value     = ft$p.value,
      n_heritable = sum(h2_cat == "Heritable"),
      n_nonher    = sum(h2_cat == "Non-heritable"),
      n_in_her    = sum(ct_vals == 1 & h2_cat == "Heritable"),
      n_in_nonher = sum(ct_vals == 1 & h2_cat == "Non-heritable")
    )
  }

  fisher_per_tissue <- list()
  for (tissue in TISSUES) {
    df <- vmr_annotated_list[[tissue]]
    ct_cols <- paste0("in_", names(atac_gr))
    for (ct_col in ct_cols) {
      res <- run_fishers(df, ct_col, tissue)
      if (!is.null(res)) fisher_per_tissue[[length(fisher_per_tissue) + 1]] <- res
    }
  }
  fisher_per_tissue_df <- bind_rows(fisher_per_tissue)
  fisher_per_tissue_df$fdr <- p.adjust(fisher_per_tissue_df$p_value, method = "fdr")

  # Pooled across tissues (deduplicated by coordinates)
  cat("Running Fisher's exact tests (pooled across tissues)...\n")
  vmr_pooled <- vmr_all |>
    distinct(seqnames, start, end, h2_category, .keep_all = TRUE)

  fisher_pooled <- list()
  ct_cols <- paste0("in_", names(atac_gr))
  for (ct_col in ct_cols) {
    res <- run_fishers(vmr_pooled, ct_col, "Pooled")
    if (!is.null(res)) fisher_pooled[[length(fisher_pooled) + 1]] <- res
  }
  fisher_pooled_df <- bind_rows(fisher_pooled)
  fisher_pooled_df$fdr <- p.adjust(fisher_pooled_df$p_value, method = "fdr")

  fisher_df <- bind_rows(fisher_per_tissue_df, fisher_pooled_df)
  fwrite(fisher_df, file.path(OUT_DIR, paste0("fishers_heritable_vs_nonheritable_", pop, ".tsv")), sep = "\t")
  cat(sprintf("Saved: fishers_heritable_vs_nonheritable_%s.tsv\n", pop))

  ## Logistic regression: in_open_chromatin ~ h2_unscaled + covariates (continuous)

  cat("\nRunning logistic regressions (continuous h²)...\n")

  run_logistic <- function(df, ct_col, tissue) {
    df_mod <- df |>
      mutate(
        y          = .data[[ct_col]],
        log_len    = log(vmr_length),
        log_snps   = log(num_snps + 1)
      )
    fit <- safe_glm(y ~ h2_unscaled + log_len + log_snps, data = df_mod)
    if (is.null(fit)) return(NULL)

    coef_df <- as.data.frame(summary(fit)$coefficients)
    coef_df$term <- rownames(coef_df)
    rownames(coef_df) <- NULL

    # Confidence intervals
    ci <- tryCatch(confint(fit), error = function(e) NULL)

    h2_row <- coef_df[coef_df$term == "h2_unscaled", ]
    if (nrow(h2_row) == 0) return(NULL)

    tibble(
      tissue    = tissue,
      cell_type = gsub("^in_", "", ct_col),
      estimate  = h2_row$Estimate,
      std_error = h2_row$`Std. Error`,
      z_value   = h2_row$`z value`,
      p_value   = h2_row$`Pr(>|z|)`,
      or        = exp(h2_row$Estimate),
      ci_lo_log = if (!is.null(ci)) ci["h2_unscaled", 1] else NA_real_,
      ci_hi_log = if (!is.null(ci)) ci["h2_unscaled", 2] else NA_real_,
      ci_lo_or  = if (!is.null(ci)) exp(ci["h2_unscaled", 1]) else NA_real_,
      ci_hi_or  = if (!is.null(ci)) exp(ci["h2_unscaled", 2]) else NA_real_,
      n_vmrs    = nrow(df_mod),
      n_in_peak = sum(df_mod$y, na.rm = TRUE)
    )
  }

  logistic_results <- list()
  for (tissue in c(TISSUES, "Pooled")) {
    df <- if (tissue == "Pooled") vmr_pooled else vmr_annotated_list[[tissue]]
    for (ct_col in paste0("in_", names(atac_gr))) {
      res <- run_logistic(df, ct_col, tissue)
      if (!is.null(res)) logistic_results[[length(logistic_results) + 1]] <- res
    }
  }
  logistic_df <- bind_rows(logistic_results)
  logistic_df$fdr <- p.adjust(logistic_df$p_value, method = "fdr")

  fwrite(logistic_df, file.path(OUT_DIR, paste0("logistic_continuous_h2_", pop, ".tsv")), sep = "\t")
  cat(sprintf("Saved: logistic_continuous_h2_%s.tsv\n", pop))

  ## Quintile analysis: open chromatin overlap rate by h² quintile

  cat("\nComputing quintile summaries...\n")

  compute_quintile <- function(df, ct_col, tissue) {
    breaks <- quantile(df$h2_unscaled,
                      probs = seq(0, 1, 1 / N_QUINTILES),
                      na.rm = TRUE) |>
      unique()
    n_bins <- length(breaks) - 1
    if (n_bins < 1) return(tibble())

    df |>
      mutate(h2_quintile = cut(h2_unscaled,
                              breaks = breaks,
                              labels = paste0("Q", seq_len(n_bins)),
                              include.lowest = TRUE)) |>
      filter(!is.na(h2_quintile)) |>
      group_by(h2_quintile) |>
      summarise(
        n          = n(),
        n_in_peak  = sum(.data[[ct_col]], na.rm = TRUE),
        frac       = n_in_peak / n,
        h2_median  = median(h2_unscaled, na.rm = TRUE),
        .groups    = "drop"
      ) |>
      mutate(
        ci_lo     = ac_ci_lo(n_in_peak, n),
        ci_hi     = ac_ci_hi(n_in_peak, n),
        cell_type = gsub("^in_", "", ct_col),
        tissue    = tissue
      )
  }

  quintile_list <- list()
  for (tissue in c(TISSUES, "Pooled")) {
    df <- if (tissue == "Pooled") vmr_pooled else vmr_annotated_list[[tissue]]
    for (ct_col in paste0("in_", names(atac_gr))) {
      res <- compute_quintile(df, ct_col, tissue)
      if (nrow(res) > 0) quintile_list[[length(quintile_list) + 1]] <- res
    }
  }
  quintile_df <- bind_rows(quintile_list)
  fwrite(quintile_df, file.path(OUT_DIR, paste0("quintile_open_chromatin_", pop, ".tsv")), sep = "\t")
  cat(sprintf("Saved: quintile_open_chromatin_%s.tsv\n", pop))

  ## Top-quintile (Q5) vs lower (Q1-Q4) Fisher's exact

  cat("\nQ5 vs Q1-Q4 Fisher's exact tests...\n")

  run_q5_fishers <- function(df, ct_col, tissue) {
    breaks <- quantile(df$h2_unscaled, probs = seq(0, 1, 0.2), na.rm = TRUE) |> unique()
    if (length(breaks) < 2) return(NULL)

    df2 <- df |>
      mutate(
        h2_quintile = cut(h2_unscaled, breaks = breaks,
                          labels = paste0("Q", seq_len(length(breaks) - 1)),
                          include.lowest = TRUE),
        is_top_q = as.integer(!is.na(h2_quintile) & h2_quintile == paste0("Q", length(breaks) - 1))
      ) |>
      filter(!is.na(h2_quintile))

    tab <- table(df2[[ct_col]], df2$is_top_q)
    if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)

    # Ensure column order: Q1-Q4 (0) vs Q5 (1)
    tab <- tab[, c("0", "1"), drop = FALSE]
    ft <- fisher.test(tab)

    tibble(
      tissue    = tissue,
      cell_type = gsub("^in_", "", ct_col),
      or        = ft$estimate,
      ci_lo     = ft$conf.int[1],
      ci_hi     = ft$conf.int[2],
      p_value   = ft$p.value,
      n_top_q   = sum(df2$is_top_q == 1),
      n_lower_q = sum(df2$is_top_q == 0)
    )
  }

  q5_list <- list()
  for (tissue in c(TISSUES, "Pooled")) {
    df <- if (tissue == "Pooled") vmr_pooled else vmr_annotated_list[[tissue]]
    for (ct_col in paste0("in_", names(atac_gr))) {
      res <- run_q5_fishers(df, ct_col, tissue)
      if (!is.null(res)) q5_list[[length(q5_list) + 1]] <- res
    }
  }
  q5_df <- bind_rows(q5_list)
  q5_df$fdr <- p.adjust(q5_df$p_value, method = "fdr")

  fwrite(q5_df, file.path(OUT_DIR, paste0("q5_vs_q1q4_fishers_", pop, ".tsv")), sep = "\t")
  cat(sprintf("Saved: q5_vs_q1q4_fishers_%s.tsv\n", pop))

  ## Summary

  cat("\n=== Fisher's exact summary (Heritable vs Non-heritable, pooled) ===\n")
  print(fisher_pooled_df |>
    dplyr::select(cell_type, or, ci_lo, ci_hi, p_value, fdr) |>
    arrange(fdr))

  cat("\n=== Logistic regression summary (continuous h², pooled) ===\n")
  print(logistic_df |>
    filter(tissue == "Pooled") |>
    dplyr::select(cell_type, or, ci_lo_or, ci_hi_or, p_value, fdr) |>
    arrange(fdr))

}

#### Reproducibility ####
cat("\nReproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
