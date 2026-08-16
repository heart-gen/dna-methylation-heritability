#### Intergenic VMR ~ Repeat Element Overlap Analysis ####
##
## Tests whether heritable VMRs — especially heritable intergenic VMRs — are
## enriched in repetitive elements / transposable elements relative to
## non-heritable VMRs.
##
## Repeat features tested:
##   Class level: LINE, SINE, LTR, DNA, Satellite, Simple_repeat, Low_complexity
##   Family level: L1, Alu, MIR, ERV1, ERVL, ERVL-MaLR, ERVK, SVA
##   Aggregate: any repeat overlap
##
## Comparisons:
##   1. All heritable vs all non-heritable VMRs (Fisher's exact)
##   2. Heritable intergenic vs non-heritable intergenic VMRs (Fisher's exact)
##   3. Continuous h² ~ repeat overlap (logistic regression + covariates)
##   4. Top h² quintile (Q5) vs Q1-Q4 within intergenic VMRs
##   5. Top h² decile vs bottom 90% within intergenic VMRs
##
## RepeatMasker source: UCSC RepeatMasker hg38 (Oct 2022) via AnnotationHub
## No liftover needed — both VMRs and RepeatMasker are in hg38.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(GenomicRanges)
  library(AnnotationHub)
})

## Configuration

TISSUES     <- c("Caudate", "DLPFC", "Hippocampus")
R2_THRESH   <- 0.3
N_QUINTILES <- 5
H2_UNIT     <- 0.1

REPEAT_CLASSES <- c("LINE", "SINE", "LTR", "DNA", "Satellite",
                    "Simple_repeat", "Low_complexity")

REPEAT_FAMILIES <- list(
  L1        = list(class = "LINE",       family = "L1"),
  Alu       = list(class = "SINE",       family = "Alu"),
  MIR       = list(class = "SINE",       family = "MIR"),
  ERV1      = list(class = "LTR",        family = "ERV1"),
  ERVL      = list(class = "LTR",        family = "ERVL"),
  ERVL_MaLR = list(class = "LTR",        family = "ERVL-MaLR"),
  ERVK      = list(class = "LTR",        family = "ERVK"),
  SVA       = list(class = "Retroposon", family = "SVA")
)

## Paths

ANNOT_DIR <- here("heritability", "elastic_net_model", "BA_only",
                  "tissue_comparison", "annotation", "_m")
ENET_BASE <- here("heritability", "elastic_net_model", "BA_only")
OUT_DIR   <- here("heritability", "elastic_net_model", "BA_only",
                  "tissue_comparison", "annotation", "repeat_elements", "_m")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## Functions

safe_glm <- function(formula, data) {
  tryCatch(
    glm(formula, data = data, family = binomial),
    error = function(e) { warning("GLM failed: ", conditionMessage(e)); NULL }
  )
}

ac_ci <- function(n_annot, n, z = qnorm(0.975)) {
  n_tilde <- n + z^2
  p_tilde <- (n_annot + z^2 / 2) / n_tilde
  lo <- pmax(0, p_tilde - z * sqrt(p_tilde * (1 - p_tilde) / n_tilde))
  hi <- pmin(1, p_tilde + z * sqrt(p_tilde * (1 - p_tilde) / n_tilde))
  list(lo = lo, hi = hi)
}

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

flag_overlaps <- function(df, annot_gr_list) {
  vmr_gr <- GRanges(
    seqnames = df$seqnames,
    ranges   = IRanges(start = df$start, end = df$end)
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

load_repeatmasker <- function() {
  cat("  Connecting to AnnotationHub...\n")
  ah   <- AnnotationHub()
  cat("  Loading AH111333 (UCSC RepeatMasker hg38, Oct 2022)...\n")
  rmsk <- ah[["AH111333"]]
  cat(sprintf("  Loaded %d repeat elements\n", length(rmsk)))
  rmsk
}

build_repeat_granges <- function(rmsk) {
  gr_list <- list()

  # Any repeat element
  cat("  Building any_repeat GRanges...\n")
  gr_list[["any_repeat"]] <- reduce(rmsk)

  # Per-class
  for (cls in REPEAT_CLASSES) {
    idx <- mcols(rmsk)$repClass == cls
    if (sum(idx) > 0) {
      gr_list[[cls]] <- reduce(rmsk[idx])
      cat(sprintf("  %s: %d elements → %d reduced regions\n",
                  cls, sum(idx), length(gr_list[[cls]])))
    }
  }

  # Per-family
  for (fam_name in names(REPEAT_FAMILIES)) {
    fam_info <- REPEAT_FAMILIES[[fam_name]]
    idx <- mcols(rmsk)$repClass  == fam_info$class &
           mcols(rmsk)$repFamily == fam_info$family
    if (sum(idx) > 0) {
      gr_list[[fam_name]] <- reduce(rmsk[idx])
      cat(sprintf("  %s (%s/%s): %d elements → %d reduced regions\n",
                  fam_name, fam_info$class, fam_info$family,
                  sum(idx), length(gr_list[[fam_name]])))
    }
  }

  gr_list
}

## Main analysis

cat("Loading RepeatMasker from AnnotationHub...\n")
rmsk      <- load_repeatmasker()
repeat_gr <- build_repeat_granges(rmsk)
cat(sprintf("\nBuilt %d repeat annotation tracks\n", length(repeat_gr)))

cat("\nRepeatMasker class composition:\n")
print(sort(table(mcols(rmsk)$repClass), decreasing = TRUE))

cat("\nLoading VMR data and computing overlaps...\n")
vmr_all_list   <- list()
vmr_inter_list <- list()

for (tissue in TISSUES) {
  cat(sprintf("\n  %s\n", tissue))
  df_all   <- load_vmrs(tissue, intergenic_only = FALSE)
  df_inter <- df_all |> filter(hg38_genes_intergenic == 1)

  cat(sprintf("    All: %d VMRs (%d Her, %d Non-her)\n",
              nrow(df_all), sum(df_all$h2_category == "Heritable"),
              sum(df_all$h2_category == "Non-heritable")))
  cat(sprintf("    Intergenic: %d VMRs (%d Her, %d Non-her)\n",
              nrow(df_inter), sum(df_inter$h2_category == "Heritable"),
              sum(df_inter$h2_category == "Non-heritable")))

  df_all   <- flag_overlaps(df_all,   repeat_gr)
  df_inter <- flag_overlaps(df_inter, repeat_gr)

  # Sanity check: any_repeat overlap rate
  for (cat_name in c("Heritable", "Non-heritable")) {
    sub_all   <- df_all[df_all$h2_category == cat_name, ]
    sub_inter <- df_inter[df_inter$h2_category == cat_name, ]
    cat(sprintf("    %s — any_repeat: all %.1f%%, intergenic %.1f%%\n",
                cat_name,
                100 * mean(sub_all$in_any_repeat),
                100 * mean(sub_inter$in_any_repeat)))
  }

  vmr_all_list[[tissue]]   <- df_all
  vmr_inter_list[[tissue]] <- df_inter
}

## Save annotated VMR table
cat("\nSaving annotated VMR table...\n")
vmr_all_out <- bind_rows(vmr_all_list)
fwrite(vmr_all_out, file.path(OUT_DIR, "vmr_repeat_overlap.tsv"), sep = "\t")
cat("Saved: vmr_repeat_overlap.tsv\n")

## Run enrichment tests

ANN_COLS_CLASS  <- paste0("in_", c("any_repeat", REPEAT_CLASSES))
ANN_COLS_FAMILY <- paste0("in_", names(REPEAT_FAMILIES))
ANN_COLS_ALL    <- c(ANN_COLS_CLASS, ANN_COLS_FAMILY)

# Fisher's exact
cat("\nRunning Fisher's exact tests...\n")
fish_list <- list()

for (tissue in TISSUES) {
  for (col in ANN_COLS_ALL) {
    res <- run_fishers(vmr_all_list[[tissue]], col, tissue, "All_VMRs")
    if (!is.null(res)) fish_list[[length(fish_list) + 1]] <- res
    res <- run_fishers(vmr_inter_list[[tissue]], col, tissue, "Intergenic_VMRs")
    if (!is.null(res)) fish_list[[length(fish_list) + 1]] <- res
  }
}

fisher_df <- bind_rows(fish_list)
fisher_df$fdr <- p.adjust(fisher_df$p_value, method = "fdr")
fwrite(fisher_df, file.path(OUT_DIR, "fishers_repeat_enrichment.tsv"), sep = "\t")
cat("Saved: fishers_repeat_enrichment.tsv\n")

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
fwrite(logistic_df, file.path(OUT_DIR, "logistic_repeat_enrichment.tsv"), sep = "\t")
cat("Saved: logistic_repeat_enrichment.tsv\n")

# Q5 and decile Fisher's exact
cat("\nRunning top-quintile and top-decile Fisher's tests (intergenic VMRs)...\n")
topbin_list <- list()

for (tissue in TISSUES) {
  df_inter <- vmr_inter_list[[tissue]]
  for (col in ANN_COLS_ALL) {
    res <- run_topbin_fishers(df_inter, col, tissue, "Intergenic_VMRs", top_frac = 0.2)
    if (!is.null(res)) topbin_list[[length(topbin_list) + 1]] <- res
    res <- run_topbin_fishers(df_inter, col, tissue, "Intergenic_VMRs", top_frac = 0.1)
    if (!is.null(res)) topbin_list[[length(topbin_list) + 1]] <- res
  }
}

topbin_df <- bind_rows(topbin_list)
topbin_df$fdr <- p.adjust(topbin_df$p_value, method = "fdr")
fwrite(topbin_df, file.path(OUT_DIR, "q5_decile_repeat_fishers.tsv"), sep = "\t")
cat("Saved: q5_decile_repeat_fishers.tsv\n")

## Summary print

cat("\n=== Fisher's enrichment: Intergenic heritable vs non-heritable (class level) ===\n")
fisher_df |>
  filter(comparison == "Intergenic_VMRs",
         annotation %in% c("any_repeat", REPEAT_CLASSES)) |>
  dplyr::select(tissue, annotation, or, ci_lo, ci_hi, p_value, fdr,
                n_in_her, n_in_nonher) |>
  arrange(annotation, tissue) |>
  print(n = Inf)

cat("\n=== Fisher's enrichment: Intergenic heritable vs non-heritable (family level) ===\n")
fisher_df |>
  filter(comparison == "Intergenic_VMRs",
         annotation %in% names(REPEAT_FAMILIES)) |>
  dplyr::select(tissue, annotation, or, ci_lo, ci_hi, p_value, fdr,
                n_in_her, n_in_nonher) |>
  arrange(annotation, tissue) |>
  print(n = Inf)

cat("\n=== Logistic regression: continuous h² ~ repeat overlap (intergenic) ===\n")
logistic_df |>
  filter(comparison == "Intergenic_VMRs") |>
  dplyr::select(tissue, annotation, or, ci_lo_or, ci_hi_or, p_value, fdr) |>
  arrange(annotation, tissue) |>
  print(n = Inf)

cat("\n=== Top-quintile / top-decile Fisher's (intergenic) ===\n")
topbin_df |>
  filter(annotation %in% c("any_repeat", "LINE", "SINE", "LTR", "L1", "Alu")) |>
  dplyr::select(tissue, comparison, annotation, or, ci_lo, ci_hi, p_value, fdr) |>
  arrange(annotation, comparison, tissue) |>
  print(n = Inf)

#### Reproducibility ####
cat("\nReproducibility information:\n")
print(Sys.time())
print(proc.time())
options(width = 120)
sessioninfo::session_info()
