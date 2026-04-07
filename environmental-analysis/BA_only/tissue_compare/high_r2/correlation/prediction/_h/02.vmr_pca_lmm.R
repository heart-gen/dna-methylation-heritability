#### VMR PCA + GLM/LMM Association Analysis (African American cohort) ####
##
## Two-stage prediction pipeline:
##   1. Load dRFE-selected VMRs per (tissue, h2_category, SDOH variable)
##   2. PCA on selected VMRs (retain PCs for cumulative variance >= 80%)
##   3a. GLM per tissue: SDOH ~ PCs + age + sex + dx + afr_ances
##   3b. LMM cross-tissue: SDOH ~ PCs + tissue + age + sex + dx + afr_ances + (1|brnum)
##
## SDOH variables are dynamically loaded from na_filter_summary.tsv (>15% NA excluded).
## Trauma variables are pooled into any_trauma_hx composite.
##
## Run from: environmental-analysis/BA_only/tissue_compare/
##           correlation/prediction/_m/

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(lme4)
  library(ordinal)
  library(broom)
  library(broom.mixed)
})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

POP_DIR      <- "BA_only"
DATA_SUFFIX  <- "AA"                 # vmr_env_assoc-{DATA_SUFFIX}.tsv.gz
TISSUES      <- c("Caudate", "DLPFC", "Hippocampus")
H2_CATS      <- c("Heritable", "Non-heritable", "Low prediction")
VAR_THRESH   <- 0.80                 # cumulative variance for PC retention
MIN_AUC      <- 0.55                 # skip tasks below chance-level dRFE score
MIN_SAMPLES  <- 40                   # minimum valid samples per task

DRFE_DIR <- here::here(
  "environmental-analysis", POP_DIR, "tissue_compare",
  "correlation", "prediction", "_m", "drfe_results"
)

COVARS <- c("age", "sex", "dx", "afr_ances")

# Variables treated as binary (binomial family) vs ordinal (clm/clmm)
BINARY_VARS  <- c("smoking", "codeine", "morphine", "cocaine", "ethanol",
                  "nicotine", "amphetamines", "any_trauma_hx", "sex", "dx")
ORDINAL_VARS <- c("education", "marital_status")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

make_file_safe <- function(s) {
  s <- tolower(trimws(s))
  s <- gsub("\\s+", "_", s)
  s <- gsub("[^\\w\\-]", "", s, perl = TRUE)
  s
}

load_pheno <- function(tissue) {
  pheno_path <- here::here("inputs", "phenotypes", "_m", "phenotypes-AA.tsv")
  ances_path <- here::here("inputs", "genetic-ancestry",
                           "structure.out_ancestry_proportion_raceDemo_compare")
  pheno <- fread(pheno_path) |>
    filter(agedeath > 17, region == tolower(tissue)) |>
    mutate(
      education = dplyr::case_when(
        education %in% c("7th","8th","Less than 7th","9th","10th","11th","12th") ~ "less_than_hs",
        education %in% c("H.S. diploma","GED") ~ "hs",
        education %in% c("1 yr college","3 yrs college","Associate's or 2 yrs college",
                         "Bachelor's","Master's","JD","PhD") ~ "more_than_hs"
      ),
      marital_status = dplyr::case_when(
        marital_status %in% c("Single") ~ "single",
        marital_status %in% c("Married") ~ "married",
        marital_status %in% c("Divorced","Separated","Widowed") ~ "previously_married"
      ),
      any_trauma_hx = dplyr::case_when(
        hx_sexual_abuse | hx_physical_abuse | hx_other_trauma | hx_military_service ~ 1L,
        is.na(hx_sexual_abuse) & is.na(hx_physical_abuse) &
          is.na(hx_other_trauma) & is.na(hx_military_service) ~ NA_integer_,
        TRUE ~ 0L
      )
    ) |>
    rename(age = agedeath, dx = primarydx) |>
    mutate_if(is.character, as.factor)

  ances <- fread(ances_path) |> filter(group == "AA")
  pheno <- left_join(pheno, ances, by = c("brnum" = "id")) |>
    rename(afr_ances = Afr)
  pheno
}

# PCA on samples x VMRs matrix (already scaled); returns list with scores and var_explained
run_pca <- function(X_scaled, var_thresh = VAR_THRESH) {
  pca <- prcomp(X_scaled, scale. = FALSE, center = FALSE)
  cum_var <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
  n_pcs <- max(2L, which(cum_var >= var_thresh)[1])
  if (is.na(n_pcs)) n_pcs <- ncol(X_scaled)
  list(
    scores      = as.data.frame(pca$x[, seq_len(n_pcs), drop = FALSE]),
    var_explained = cum_var[n_pcs],
    n_pcs       = n_pcs
  )
}

# Impute a matrix column-wise with column medians
impute_median <- function(X) {
  apply(X, 2, function(col) {
    med <- median(col, na.rm = TRUE)
    ifelse(is.na(col), med, col)
  })
}

# Fit GLM or CLM and return tidy coefficient table for PCs only
fit_model_tissue <- function(model_df, env_var, pc_cols) {
  formula_str <- paste0(env_var, " ~ ", paste(pc_cols, collapse = " + "),
                        " + age + sex + dx + afr_ances")
  tryCatch({
    if (env_var %in% BINARY_VARS) {
      model_df[[env_var]] <- as.integer(model_df[[env_var]])
      fit <- glm(as.formula(formula_str), data = model_df, family = binomial)
      tidy(fit) |> filter(grepl("^PC", term))
    } else {
      model_df[[env_var]] <- factor(model_df[[env_var]],
                                    levels = c("less_than_hs","hs","more_than_hs",
                                               "single","married","previously_married"),
                                    ordered = TRUE)
      fit <- clm(as.formula(formula_str), data = model_df)
      tidy(fit) |> filter(grepl("^PC", term))
    }
  }, error = function(e) {
    warning(sprintf("  GLM failed for %s: %s", env_var, conditionMessage(e)))
    NULL
  })
}

# Fit GLMM or CLMM for cross-tissue LMM and return tidy PC rows
fit_model_crosst <- function(model_df, env_var, pc_cols) {
  formula_str <- paste0(
    env_var, " ~ ", paste(pc_cols, collapse = " + "),
    " + tissue + age + sex + dx + afr_ances + (1|brnum)"
  )
  tryCatch({
    if (env_var %in% BINARY_VARS) {
      model_df[[env_var]] <- as.integer(model_df[[env_var]])
      fit <- glmer(as.formula(formula_str), data = model_df, family = binomial,
                   control = glmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 1e5)))
      tidy(fit, effects = "fixed") |> filter(grepl("^PC", term))
    } else {
      model_df[[env_var]] <- factor(model_df[[env_var]],
                                    levels = c("less_than_hs","hs","more_than_hs",
                                               "single","married","previously_married"),
                                    ordered = TRUE)
      fit <- clmm(as.formula(formula_str), data = model_df)
      tidy(fit) |> filter(grepl("^PC", term))
    }
  }, error = function(e) {
    warning(sprintf("  LMM failed for %s: %s", env_var, conditionMessage(e)))
    NULL
  })
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

## Load NA filter results to get included SDOH variables
na_filter_path <- file.path(DRFE_DIR, "na_filter_summary.tsv")
if (!file.exists(na_filter_path)) {
  stop("na_filter_summary.tsv not found. Run 00.na_filter_analysis.py first.")
}
na_filter  <- read.delim(na_filter_path)
ENV_VARS   <- na_filter$variable[!na_filter$excluded]
cat("SDOH variables included:", paste(ENV_VARS, collapse = ", "), "\n\n")

## Output directories
out_dir <- file.path(DRFE_DIR)

all_results_tissue <- list()
all_results_crosst <- list()

for (tissue in TISSUES) {
  cat("===", tissue, "===\n")

  drfe_summary_path <- file.path(DRFE_DIR, paste0("drfe_summary_", tissue, ".tsv"))
  if (!file.exists(drfe_summary_path)) {
    warning("  No dRFE summary for ", tissue, " — skipping")
    next
  }
  drfe_sum <- read.delim(drfe_summary_path)

  data_path <- here::here(
    "environmental-analysis", POP_DIR, tolower(tissue),
    "correlation", "_m", paste0("vmr_env_assoc-", DATA_SUFFIX, ".tsv.gz")
  )
  if (!file.exists(data_path)) {
    warning("  Data file not found: ", data_path)
    next
  }
  cat("  Loading methylation data...\n")
  df <- fread(data_path)

  pheno_tissue <- load_pheno(tissue)

  for (h2_cat in H2_CATS) {
    h2_safe <- make_file_safe(h2_cat)
    df_sub  <- df[df$h2_category == h2_cat, ]

    X_wide <- df_sub |>
      dplyr::select(brnum, feature_id, meth) |>
      pivot_wider(names_from = feature_id, values_from = meth)

    if (nrow(X_wide) < MIN_SAMPLES || ncol(X_wide) < 3) next

    for (env_var in ENV_VARS) {
      task_row <- drfe_sum[drfe_sum$h2_category == h2_cat &
                             drfe_sum$env_var == env_var, ]
      if (nrow(task_row) == 0) next
      if (task_row$best_score[1] < MIN_AUC) next

      best_n <- task_row$best_n_features[1]
      if (is.na(best_n) || best_n < 1) next

      feat_path <- file.path(DRFE_DIR,
                             paste0(tissue, "_", h2_safe, "_", env_var),
                             "ranked_features.tsv")
      if (!file.exists(feat_path)) next

      ranked <- read.delim(feat_path)
      selected_vmrs <- head(ranked$feature, best_n)

      vmr_cols <- intersect(selected_vmrs, colnames(X_wide))
      if (length(vmr_cols) < 2) next

      X_mat <- as.matrix(X_wide[, vmr_cols])
      rownames(X_mat) <- X_wide$brnum

      X_mat    <- impute_median(X_mat)
      X_scaled <- scale(X_mat)

      ok_cols <- apply(X_scaled, 2, function(x) var(x, na.rm = TRUE) > 0)
      X_scaled <- X_scaled[, ok_cols, drop = FALSE]
      if (ncol(X_scaled) < 2) next

      pca_res <- run_pca(X_scaled)
      pcs_df  <- pca_res$scores
      pcs_df$brnum <- X_wide$brnum
      pc_cols <- paste0("PC", seq_len(pca_res$n_pcs))

      model_df <- pcs_df |>
        left_join(pheno_tissue, by = "brnum") |>
        drop_na(all_of(c(env_var, COVARS)))

      if (nrow(model_df) < MIN_SAMPLES) next
      if (length(unique(model_df[[env_var]])) < 2) next

      cat(sprintf("  %s / %s / %s: n=%d, n_vmrs=%d, n_pcs=%d, var_exp=%.1f%%\n",
                  tissue, h2_cat, env_var, nrow(model_df),
                  length(vmr_cols), pca_res$n_pcs,
                  pca_res$var_explained * 100))

      coef_df <- fit_model_tissue(model_df, env_var, pc_cols)
      if (!is.null(coef_df) && nrow(coef_df) > 0) {
        coef_df$fdr <- p.adjust(coef_df$p.value, method = "fdr")
        coef_df <- coef_df |>
          mutate(tissue = tissue, h2_category = h2_cat,
                 env_var = env_var, n_samples = nrow(model_df),
                 n_pcs = pca_res$n_pcs, var_explained = pca_res$var_explained)
        all_results_tissue[[length(all_results_tissue) + 1]] <- coef_df
      }
    }
  }
}

## ----- Cross-tissue LMM -----
cat("\n=== Cross-tissue LMM ===\n")

for (h2_cat in H2_CATS) {
  h2_safe <- make_file_safe(h2_cat)

  for (env_var in ENV_VARS) {
    tissue_vmrs <- list()
    tissue_dfs  <- list()

    for (tissue in TISSUES) {
      drfe_sum_path <- file.path(DRFE_DIR, paste0("drfe_summary_", tissue, ".tsv"))
      if (!file.exists(drfe_sum_path)) next
      drfe_sum  <- read.delim(drfe_sum_path)
      task_row  <- drfe_sum[drfe_sum$h2_category == h2_cat &
                              drfe_sum$env_var == env_var, ]
      if (nrow(task_row) == 0 || task_row$best_score[1] < MIN_AUC) next

      feat_path <- file.path(DRFE_DIR,
                             paste0(tissue, "_", h2_safe, "_", env_var),
                             "ranked_features.tsv")
      if (!file.exists(feat_path)) next

      ranked <- read.delim(feat_path)
      tissue_vmrs[[tissue]] <- head(ranked$feature, task_row$best_n_features[1])

      data_path <- here::here(
        "environmental-analysis", POP_DIR, tolower(tissue),
        "correlation", "_m", paste0("vmr_env_assoc-", DATA_SUFFIX, ".tsv.gz")
      )
      if (!file.exists(data_path)) next
      df <- fread(data_path)
      df_sub <- df[df$h2_category == h2_cat, ]

      X_wide <- df_sub |>
        dplyr::select(brnum, feature_id, meth) |>
        pivot_wider(names_from = feature_id, values_from = meth)

      tissue_dfs[[tissue]] <- X_wide
    }

    if (length(tissue_vmrs) < 2) next

    union_vmrs <- unique(unlist(tissue_vmrs))
    combined_list <- lapply(names(tissue_dfs), function(tis) {
      X_w <- tissue_dfs[[tis]]
      vmr_cols <- intersect(union_vmrs, colnames(X_w))
      if (length(vmr_cols) < 2) return(NULL)

      X_mat <- as.matrix(X_w[, vmr_cols])
      rownames(X_mat) <- X_w$brnum

      miss_vmrs <- setdiff(union_vmrs, vmr_cols)
      if (length(miss_vmrs) > 0) {
        na_mat <- matrix(NA_real_, nrow = nrow(X_mat), ncol = length(miss_vmrs),
                         dimnames = list(rownames(X_mat), miss_vmrs))
        X_mat <- cbind(X_mat, na_mat)
      }
      X_mat <- X_mat[, union_vmrs, drop = FALSE]

      X_mat    <- impute_median(X_mat)
      X_scaled <- scale(X_mat)
      df_out   <- as.data.frame(X_scaled)
      df_out$brnum  <- X_w$brnum
      df_out$tissue <- tis
      df_out
    })
    combined_list <- Filter(Negate(is.null), combined_list)
    if (length(combined_list) < 2) next

    combined <- bind_rows(combined_list)

    vmr_subset <- setdiff(colnames(combined), c("brnum", "tissue"))
    ok_cols    <- sapply(combined[, vmr_subset], function(x) var(x, na.rm = TRUE) > 0)
    vmr_subset <- vmr_subset[ok_cols]
    if (length(vmr_subset) < 2) next

    X_combined <- as.matrix(combined[, vmr_subset])

    pca_res <- run_pca(X_combined)
    pcs_df  <- pca_res$scores
    pcs_df$brnum  <- combined$brnum
    pcs_df$tissue <- combined$tissue
    pc_cols <- paste0("PC", seq_len(pca_res$n_pcs))

    pheno_all <- lapply(TISSUES, function(tis) {
      p <- load_pheno(tis)
      p$tissue <- tis
      p
    })
    pheno_all <- bind_rows(pheno_all)

    model_df <- pcs_df |>
      left_join(pheno_all, by = c("brnum", "tissue")) |>
      drop_na(all_of(c(env_var, COVARS)))

    if (nrow(model_df) < MIN_SAMPLES) next
    if (length(unique(model_df[[env_var]])) < 2) next

    cat(sprintf("  %s / %s: n_obs=%d, n_tissues=%d, n_pcs=%d, var_exp=%.1f%%\n",
                h2_cat, env_var, nrow(model_df),
                length(unique(model_df$tissue)), pca_res$n_pcs,
                pca_res$var_explained * 100))

    coef_df <- fit_model_crosst(model_df, env_var, pc_cols)
    if (!is.null(coef_df) && nrow(coef_df) > 0) {
      coef_df$fdr <- p.adjust(coef_df$p.value, method = "fdr")
      coef_df <- coef_df |>
        mutate(h2_category = h2_cat, env_var = env_var,
               n_obs = nrow(model_df),
               n_tissues = length(unique(model_df$tissue)),
               n_pcs = pca_res$n_pcs,
               var_explained = pca_res$var_explained)
      all_results_crosst[[length(all_results_crosst) + 1]] <- coef_df
    }
  }
}

## ----- Save outputs -----
if (length(all_results_tissue) > 0) {
  res_tissue <- bind_rows(all_results_tissue)
  fwrite(res_tissue,
         file.path(out_dir, "pca_glm_tissue_results.tsv"),
         sep = "\t", quote = FALSE)
  cat("\nPer-tissue GLM results saved:", nrow(res_tissue), "rows\n")
}

if (length(all_results_crosst) > 0) {
  res_crosst <- bind_rows(all_results_crosst)
  fwrite(res_crosst,
         file.path(out_dir, "pca_lmm_crosst_results.tsv"),
         sep = "\t", quote = FALSE)
  cat("Cross-tissue LMM results saved:", nrow(res_crosst), "rows\n")
}

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
