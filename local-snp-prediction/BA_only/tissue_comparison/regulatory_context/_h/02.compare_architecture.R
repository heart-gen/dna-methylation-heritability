#### Architecture-adjusted comparison of VMR regulatory-context associations ####
#
# Args: cohort, tissue, modality, population, window,
#       expression_link_layers (optional; expression only: both | abc | nearest_gene).
# PSI uses window_*kb only. Expression defaults to both abc and nearest_gene_window_*kb.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(broom)
})

source(here("heritability", "elastic_net_model", "BA_only",
            "tissue_comparison", "regulatory_context", "_h",
            "00.regulatory_context_utils.R"))

args       <- commandArgs(trailingOnly = TRUE)
cohort     <- ifelse(length(args) >= 1, args[[1]], "BA_only")
tissue     <- ifelse(length(args) >= 2, tolower(args[[2]]), "dlpfc")
modality   <- ifelse(length(args) >= 3, tolower(args[[3]]), "expression")
population <- ifelse(length(args) >= 4, toupper(args[[4]]), "AA")
window     <- ifelse(length(args) >= 5, as.integer(args[[5]]), 250000L)
expression_layers <- ifelse(length(args) >= 6, args[[6]], "both")

run_tags <- resolve_regulatory_run_tags(modality, window, expression_layers)

for (run_tag in run_tags) {
  in_dir  <- here("heritability", "elastic_net_model", cohort,
                  "tissue_comparison", "regulatory_context", "_m",
                  tissue, population, modality, run_tag)
  out_dir <- in_dir

  summary_fn <- file.path(in_dir, "vmr_association_summary.tsv")
  links_fn   <- file.path(in_dir, "tested_vmr_feature_links.tsv.gz")
  if (!file.exists(summary_fn) || !file.exists(links_fn)) {
    message2(
      "Skipping layer run_tag=%s — missing inputs under %s",
      run_tag, in_dir
    )
    next
  }

  assoc_summary <- fread(summary_fn)
  links <- fread(links_fn)
  pheno <- load_pheno(cohort, tissue)
  meth  <- load_methylation(cohort, tissue, pheno)

  meth_var <- meth |> group_by(vmr_id) |>
    summarise(methylation_variance = var(meth, na.rm = TRUE), .groups = "drop")

  distance_df <- links |> group_by(vmr_id) |>
    summarise(min_distance = min(distance, na.rm = TRUE), .groups = "drop") |>
    mutate(min_distance = ifelse(is.infinite(min_distance), NA_real_, min_distance))

  arch <- load_architecture_covariates(cohort, tissue, population) |>
    dplyr::select(vmr_id, vmr_length, num_snps, chromatin_state,
                  cpg_density_proxy, h2_category, h2_unscaled, r_squared_cv)

  df <- assoc_summary |>
    dplyr::select(-any_of(c("h2_category", "h2_unscaled"))) |>
    left_join(arch, by = "vmr_id") |>
    left_join(meth_var, by = "vmr_id") |>
    left_join(distance_df, by = "vmr_id") |>
    mutate(
      h2_category = factor(h2_category, levels = H2_GROUP_LEVELS),
      any_sig_fdr_05 = as.integer(n_sig_fdr_05 > 0),
      any_sig_fdr_10 = as.integer(n_sig_fdr_10 > 0),
      h2_scaled = h2_unscaled / 0.1,
      log_vmr_length = log(pmax(vmr_length, 1)),
      log_num_snps = log(pmax(num_snps, 0) + 1),
      log_min_distance = log(pmax(min_distance, 0) + 1),
      log_max_abs_beta = log(pmax(max_abs_beta, 0) + 1e-8)
    )

  safe_fwrite(df, file.path(out_dir, "architecture_model_input.tsv"), sep = "\t")

  complete_vars <- c("h2_category", "log_vmr_length", "log_num_snps",
                     "methylation_variance", "log_min_distance",
                     "chromatin_state")
  df_cat <- df |>
    filter(!is.na(h2_category), complete.cases(across(all_of(complete_vars))))

  safe_glm_binom <- function(formula, data) {
    tryCatch(glm(formula, data = data, family = binomial),
             error = function(e) {
               warning("GLM failed: ", conditionMessage(e))
               NULL
             })
  }

  safe_lm_gauss <- function(formula, data) {
    tryCatch(lm(formula, data = data),
             error = function(e) {
               warning("LM failed: ", conditionMessage(e))
               NULL
             })
  }

  model_rows <- list()
  if (nrow(df_cat) > 20 && length(unique(df_cat$any_sig_fdr_10)) > 1) {
    fit_burden <- safe_glm_binom(
      any_sig_fdr_10 ~ h2_category + log_vmr_length + log_num_snps +
        methylation_variance + log_min_distance + chromatin_state,
      df_cat
    )
    if (!is.null(fit_burden)) {
      model_rows$burden_category <- tidy(fit_burden, conf.int = TRUE,
                                         exponentiate = TRUE) |>
        mutate(model = "burden_fdr10_category", estimate_scale = "odds_ratio")
    }
  }

  df_strength <- df_cat |> filter(is.finite(log_max_abs_beta))
  if (nrow(df_strength) > 20) {
    fit_strength <- safe_lm_gauss(
      log_max_abs_beta ~ h2_category + log_vmr_length + log_num_snps +
        methylation_variance + log_min_distance + chromatin_state,
      df_strength
    )
    if (!is.null(fit_strength)) {
      model_rows$strength_category <- tidy(fit_strength, conf.int = TRUE) |>
        mutate(model = "strength_category", estimate_scale = "log_abs_beta")
    }
  }

  df_cont <- df |>
    filter(r_squared_cv > 0.3, !is.na(h2_unscaled),
           complete.cases(across(all_of(setdiff(complete_vars, "h2_category")))))
  if (nrow(df_cont) > 20 && length(unique(df_cont$any_sig_fdr_10)) > 1) {
    fit_cont <- safe_glm_binom(
      any_sig_fdr_10 ~ h2_scaled + log_vmr_length + log_num_snps +
        methylation_variance + log_min_distance + chromatin_state,
      df_cont
    )
    if (!is.null(fit_cont)) {
      model_rows$burden_continuous <- tidy(fit_cont, conf.int = TRUE,
                                           exponentiate = TRUE) |>
        mutate(model = "burden_fdr10_continuous_h2",
               estimate_scale = "odds_ratio_per_0.1_h2")
    }
  }

  model_results <- bind_rows(model_rows)
  if (nrow(model_results) > 0) {
    model_results <- model_results |>
      mutate(
        cohort = cohort, tissue = tissue, modality = modality,
        population = population, run_tag = run_tag
      )
  }
  safe_fwrite(model_results, file.path(out_dir, "architecture_adjusted_models.tsv"),
              sep = "\t")

  make_bins <- function(x, n_bins, prefix) {
    br <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1),
                          na.rm = TRUE))
    if (length(br) < 2) return(factor(rep(NA_character_, length(x))))
    cut(x, breaks = br, include.lowest = TRUE,
        labels = paste0(prefix, seq_len(length(br) - 1)))
  }

  bin_summary <- df_cont |>
    mutate(
      h2_tertile = make_bins(h2_unscaled, 3, "T"),
      h2_quintile = make_bins(h2_unscaled, 5, "Q")
    ) |>
    pivot_longer(c(h2_tertile, h2_quintile),
                 names_to = "bin_type", values_to = "h2_bin") |>
    filter(!is.na(h2_bin)) |>
    group_by(cohort = cohort, tissue = tissue, population = population,
             modality = modality, run_tag = run_tag, bin_type, h2_bin) |>
    summarise(
      n_vmrs = n(),
      n_any_fdr10 = sum(any_sig_fdr_10, na.rm = TRUE),
      prop_any_fdr10 = mean(any_sig_fdr_10, na.rm = TRUE),
      median_max_abs_beta = median(max_abs_beta, na.rm = TRUE),
      median_pairs_tested = median(n_pairs_tested, na.rm = TRUE),
      .groups = "drop"
    )
  safe_fwrite(bin_summary, file.path(out_dir, "h2_bin_association_summary.tsv"),
              sep = "\t")

  group_summary <- df |>
    filter(!is.na(h2_category)) |>
    group_by(cohort = cohort, tissue = tissue, population = population,
             modality = modality, run_tag = run_tag, h2_category) |>
    summarise(
      n_vmrs = n(),
      n_any_fdr05 = sum(any_sig_fdr_05, na.rm = TRUE),
      n_any_fdr10 = sum(any_sig_fdr_10, na.rm = TRUE),
      prop_any_fdr05 = mean(any_sig_fdr_05, na.rm = TRUE),
      prop_any_fdr10 = mean(any_sig_fdr_10, na.rm = TRUE),
      median_max_abs_beta = median(max_abs_beta, na.rm = TRUE),
      median_pairs_tested = median(n_pairs_tested, na.rm = TRUE),
      .groups = "drop"
    )
  safe_fwrite(group_summary, file.path(out_dir, "h2_group_association_summary.tsv"),
              sep = "\t")

  message2("Saved architecture comparison outputs to %s [run_tag=%s]",
           out_dir, run_tag)
}

#### Reproducibility ####
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
