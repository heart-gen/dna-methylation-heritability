#### Continuous h² ~ Genomic Annotation GLM Analysis ####
##
## Sensitivity / refinement analysis: test whether SNP-based heritability (h²)
## estimated from boosting elastic net is monotonically associated with genomic
## annotation membership across VMRs.
##
## This is NOT a replacement for the manuscript's 3-class framework
## (Heritable / Non-heritable / Low prediction). It asks: within the space of
## high-confidence VMRs (r² > 0.3), does annotation composition shift
## continuously with h²?
##
## Primary analysis: VMRs with r_squared_cv > 0.3
## Sensitivity analysis: all VMRs with valid h² (including low-prediction)
##
## Models per annotation × tissue:
##   1. Linear logistic: annotation ~ h2_scaled + log(vmr_length) + log(num_snps+1)
##   2. Spline logistic: annotation ~ ns(h2_scaled, df=3) + confounders (LRT vs null)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(splines)
  library(broom)
})

## Configuration

TISSUES    <- c("Caudate", "DLPFC", "Hippocampus")
POPULATIONS <- c("AA", "EA") # using matched VMRs, population-specific h2 values
R2_THRESH  <- 0.3          # exclude low-prediction VMRs in primary analysis
N_QUINTILES <- 5
H2_UNIT    <- 0.1          # OR reported per this unit increase in h2_unscaled

ANNOT_DIR <- here::here(
  "heritability", "elastic_net_model", "all_individuals",
  "tissue_comparison", "annotation", "_m"
)
ENET_BASE <- here::here(
  "heritability", "elastic_net_model", "all_individuals"
)
OUT_DIR <- here::here(
  "heritability", "elastic_net_model", "all_individuals",
  "tissue_comparison", "annotation", "h2_continuous", "_m"
)

# Annotation columns to test → display labels
ANNOT_COLS <- c(
  "hg38_genes_promoters"  = "Promoter",
  "hg38_enhancers_fantom" = "Enhancer",
  "hg38_genes_1to5kb"     = "1\u20135 kb upstream",
  "hg38_genes_intergenic" = "Intergenic"
)

## Helpers

load_data <- function(tissue, pop) {
  tissue_lower <- tolower(tissue)

  # Wide annotation: one row per VMR, binary 0/1 annotation columns
  annot <- fread(file.path(ANNOT_DIR,
    paste0(tissue_lower, "_vmr_annotations_hg38_wide.tsv")))

  # Elastic-net summary: h2_unscaled, r_squared_cv, num_snps
  enet <- fread(file.path(ENET_BASE, tissue_lower, "_m",
    paste0(tissue_lower, "_summary_elastic-net_matched_r2_0.3.tsv"))) |>
    rename(seqnames = "chrom") |>
    dplyr::select(seqnames, start, end,
                  h2_unscaled = paste0("h2_unscaled_", pop),
                  r_squared_cv = paste0("r_squared_cv_", pop),
                  num_snps = paste0("num_snps_", pop)) |>
    mutate(.enet_matched = TRUE)

  # Join on genomic coordinates
  df <- annot |>
    left_join(enet, by = c("seqnames", "start", "end")) |>
    mutate(
      vmr_length = end - start,
      h2_scaled  = h2_unscaled / H2_UNIT,   # OR per 0.1 h2 unit
      tissue     = tissue
    )

  n_unmatched <- sum(is.na(df$.enet_matched))
  if (n_unmatched > 0) {
    warning(sprintf(
      "%s: %d/%d annotation rows lacked elastic-net matches after coordinate join",
      tissue, n_unmatched, nrow(df)
    ))
  }

  df |>
    dplyr::select(-.enet_matched)
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

# Agresti-Coull 95% CI bounds (vectorized, safe for use inside mutate)
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

# Linear logistic regression: OR per H2_UNIT increase
run_glm_linear <- function(df, annot_col, tissue) {
  fit <- safe_glm(
    as.formula(paste0(
      annot_col,
      " ~ h2_scaled + log(vmr_length) + log(num_snps + 1)"
    )),
    data = df
  )
  if (is.null(fit)) return(NULL)

  tidy(fit, conf.int = TRUE, exponentiate = TRUE) |>
    filter(term == "h2_scaled") |>
    mutate(
      annotation  = ANNOT_COLS[[annot_col]],
      annot_col   = annot_col,
      tissue      = tissue,
      n_vmrs      = nrow(df),
      n_annotated = sum(df[[annot_col]], na.rm = TRUE)
    )
}

# Spline logistic regression: LRT of spline vs covariate-only null
run_glm_spline_lrt <- function(df, annot_col, tissue) {
  fit_spline <- safe_glm(
    as.formula(paste0(
      annot_col,
      " ~ ns(h2_scaled, df = 3) + log(vmr_length) + log(num_snps + 1)"
    )),
    data = df
  )
  fit_null <- safe_glm(
    as.formula(paste0(
      annot_col,
      " ~ log(vmr_length) + log(num_snps + 1)"
    )),
    data = df
  )
  if (is.null(fit_spline) || is.null(fit_null)) return(NULL)

  lrt <- anova(fit_null, fit_spline, test = "LRT")
  tibble(
    annotation  = ANNOT_COLS[[annot_col]],
    annot_col   = annot_col,
    tissue      = tissue,
    df_diff     = lrt$Df[2],
    deviance    = lrt$Deviance[2],
    p_lrt       = lrt$`Pr(>Chi)`[2],
    n_vmrs      = nrow(df),
    n_annotated = sum(df[[annot_col]], na.rm = TRUE)
  )
}

# Spline predictions over a h2 grid at median confounder values
compute_spline_predictions <- function(df, annot_col, tissue) {
  fit_spline <- safe_glm(
    as.formula(paste0(
      annot_col,
      " ~ ns(h2_scaled, df = 3) + log(vmr_length) + log(num_snps + 1)"
    )),
    data = df
  )
  if (is.null(fit_spline)) return(NULL)

  h2_grid <- seq(
    quantile(df$h2_unscaled, 0.02, na.rm = TRUE),
    quantile(df$h2_unscaled, 0.98, na.rm = TRUE),
    length.out = 200
  )
  newdata <- data.frame(
    h2_unscaled = h2_grid,
    h2_scaled   = h2_grid / H2_UNIT,
    vmr_length  = median(df$vmr_length, na.rm = TRUE),
    num_snps    = median(df$num_snps,   na.rm = TRUE)
  )
  pred <- predict(fit_spline, newdata = newdata, type = "link", se.fit = TRUE)
  tibble(
    h2_unscaled = h2_grid,
    prob        = plogis(pred$fit),
    ci_lo       = plogis(pred$fit - 1.96 * pred$se.fit),
    ci_hi       = plogis(pred$fit + 1.96 * pred$se.fit),
    annotation  = ANNOT_COLS[[annot_col]],
    annot_col   = annot_col,
    tissue      = tissue
  )
}

# Annotation fraction per h2 quintile with Agresti-Coull CI
compute_quintile_summary <- function(df, annot_col, tissue) {
  breaks <- quantile(df$h2_unscaled,
                     probs = seq(0, 1, 1 / N_QUINTILES),
                     na.rm = TRUE) |>
    unique()
  n_bins <- length(breaks) - 1

  if (n_bins < 1) {
    warning(sprintf(
      "%s / %s: insufficient unique h2 values to compute quantile bins",
      tissue, ANNOT_COLS[[annot_col]]
    ))
    return(tibble())
  }

  if (n_bins < N_QUINTILES) {
    warning(sprintf(
      "%s / %s: reduced h2 bins from %d to %d because quantile breaks were not unique",
      tissue, ANNOT_COLS[[annot_col]], N_QUINTILES, n_bins
    ))
  }

  df |>
    mutate(h2_quintile = cut(h2_unscaled,
                             breaks = breaks,
                             labels = paste0("Q", seq_len(n_bins)),
                             include.lowest = TRUE)) |>
    filter(!is.na(h2_quintile)) |>
    group_by(h2_quintile) |>
    summarise(
      n          = n(),
      n_annotated = sum(.data[[annot_col]], na.rm = TRUE),
      frac       = n_annotated / n,
      h2_median  = median(h2_unscaled, na.rm = TRUE),
      h2_lo      = min(h2_unscaled,    na.rm = TRUE),
      h2_hi      = max(h2_unscaled,    na.rm = TRUE),
      .groups    = "drop"
    ) |>
    mutate(
      ci_lo       = ac_ci_lo(n_annotated, n),
      ci_hi       = ac_ci_hi(n_annotated, n),
      annotation  = ANNOT_COLS[[annot_col]],
      annot_col   = annot_col,
      tissue      = tissue
    )
}

## Main

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (pop in POPULATIONS){

  all_glm      <- list()
  all_lrt      <- list()
  all_spline   <- list()
  all_quint    <- list()
  all_glm_sa   <- list()
  all_quint_sa <- list()

  for (tissue in TISSUES) {
    cat("===", tissue, "===\n")

    df_all <- load_data(tissue, pop)

    # Primary: high-confidence VMRs only
    df <- df_all |>
      filter(
        !is.na(r_squared_cv), r_squared_cv > R2_THRESH,
        !is.na(h2_unscaled),  is.finite(h2_unscaled),
        vmr_length > 0, !is.na(num_snps), num_snps > 0
      )
    cat(sprintf("  Primary (r2 > %.1f): %d VMRs\n", R2_THRESH, nrow(df)))

    # Sensitivity: all VMRs with valid h2 (include low-prediction)
    df_sa <- df_all |>
      filter(
        !is.na(h2_unscaled),  is.finite(h2_unscaled),
        !is.na(r_squared_cv),
        vmr_length > 0, !is.na(num_snps), num_snps > 0
      )
    cat(sprintf("  Sensitivity (all valid h2): %d VMRs\n", nrow(df_sa)))

    for (ac in names(ANNOT_COLS)) {
      if (!ac %in% colnames(df)) {
        warning("  Column not found: ", ac, " — skipping")
        next
      }

      # Primary analysis
      glm_res  <- run_glm_linear(df, ac, tissue)
      lrt_res  <- run_glm_spline_lrt(df, ac, tissue)
      spl_res  <- compute_spline_predictions(df, ac, tissue)
      qt_res   <- compute_quintile_summary(df, ac, tissue)

      if (!is.null(glm_res)) all_glm[[length(all_glm) + 1]]     <- glm_res
      if (!is.null(lrt_res)) all_lrt[[length(all_lrt) + 1]]     <- lrt_res
      if (!is.null(spl_res)) all_spline[[length(all_spline) + 1]] <- spl_res
      all_quint[[length(all_quint) + 1]] <- qt_res

      # Sensitivity analysis
      glm_sa <- run_glm_linear(df_sa, ac, tissue)
      qt_sa  <- compute_quintile_summary(df_sa, ac, tissue)

      if (!is.null(glm_sa)) all_glm_sa[[length(all_glm_sa) + 1]] <- glm_sa
      all_quint_sa[[length(all_quint_sa) + 1]] <- qt_sa
    }
  }

  # FDR correction across all annotation × tissue tests
  fdr_adjust <- function(lst, p_col = "p.value") {
    df <- bind_rows(lst)
    df$fdr <- p.adjust(df[[p_col]], method = "fdr")
    df
  }

  res_glm      <- fdr_adjust(all_glm)
  res_lrt      <- fdr_adjust(all_lrt, p_col = "p_lrt")
  res_spline   <- bind_rows(all_spline)
  res_quint    <- bind_rows(all_quint)
  res_glm_sa   <- fdr_adjust(all_glm_sa)
  res_quint_sa <- bind_rows(all_quint_sa)

  fwrite(res_glm,      file.path(OUT_DIR, paste0("glm_linear_results_", pop, ".tsv")),          sep = "\t")
  fwrite(res_lrt,      file.path(OUT_DIR, paste0("spline_lrt_results_", pop, ".tsv")),          sep = "\t")
  fwrite(res_spline,   file.path(OUT_DIR, paste0("spline_predictions_", pop, ".tsv")),          sep = "\t")
  fwrite(res_quint,    file.path(OUT_DIR, paste0("quintile_summary_", pop, ".tsv")),            sep = "\t")
  fwrite(res_glm_sa,   file.path(OUT_DIR, paste0("glm_linear_sensitivity_", pop, ".tsv")),     sep = "\t")
  fwrite(res_quint_sa, file.path(OUT_DIR, paste0("quintile_summary_sensitivity_", pop, ".tsv")), sep = "\t")

  cat("\nResults written to:", OUT_DIR, "\n")
  cat("Linear GLM summary:\n")
  print(res_glm |>
    dplyr::select(tissue, annotation, estimate, conf.low, conf.high, p.value, fdr) |>
    arrange(fdr))
}


#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
