#### Environment/clinical convergence among DNAm-linked transcriptional features ####
#
# Args: cohort, tissue, modality, population, window,
#       max_features (optional smoke-test limit),
#       expression_link_layers (optional; expression only: both | abc | nearest_gene),
#       vmr_set (optional: shared | AA_only | EA_only).

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(limma)
})

source(here("heritability", "elastic_net_model", "all_individuals",
            "tissue_comparison", "regulatory_context", "_h",
            "00.regulatory_context_utils.R"))
source(here("environmental-analysis", "BA_only", "tissue_compare",
            "_h", "discovery_exposures.R"))

args <- commandArgs(trailingOnly = TRUE)
cohort <- ifelse(length(args) >= 1, args[[1]], "all_individuals")
tissue <- ifelse(length(args) >= 2, tolower(args[[2]]), "dlpfc")
modality <- ifelse(length(args) >= 3, tolower(args[[3]]), "expression")
population <- ifelse(length(args) >= 4, toupper(args[[4]]), "AA")
window <- ifelse(length(args) >= 5, as.integer(args[[5]]), 250000L)
max_features <- if (length(args) >= 6 &&
                    !toupper(args[[6]]) %in% c("", "NA", "NULL", "NONE")) {
  as.integer(args[[6]])
} else {
  NA_integer_
}
expression_layers <- ifelse(length(args) >= 7, args[[7]], "both")
vmr_set <- ifelse(length(args) >= 8, args[[8]], "shared")
vmr_set <- validate_vmr_set(cohort, population, vmr_set)

if (should_skip_shared_duplicate_population(population, vmr_set)) {
  message2(
    paste0(
      "Skipping duplicate shared-VMR environment convergence for population=%s ",
      "(uses canonical association inputs from population=%s)."
    ),
    population, SHARED_VMR_CANONICAL_POPULATION
  )
  quit(save = "no", status = 0)
}

run_tags <- resolve_regulatory_run_tags(modality, window, expression_layers)

for (run_tag in run_tags) {
  assoc_pop <- regctx_assoc_source_population(population, vmr_set)
  in_dir <- regctx_output_dir(
    cohort, tissue, assoc_pop, modality, run_tag, vmr_set
  )
  out_dir <- regctx_output_dir(
    cohort, tissue, population, modality, run_tag, vmr_set
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  assoc_fn <- file.path(in_dir, "vmr_feature_associations.tsv.gz")
  if (!file.exists(assoc_fn)) {
    message2(
      "Skipping environment convergence for run_tag=%s — missing %s",
      run_tag, assoc_fn
    )
    next
  }

  assoc <- fread(assoc_fn)
  assoc <- ensure_vmr_set_column(assoc, vmr_set)
  tested_features <- unique(assoc$feature_id)
  if (!is.na(max_features) && length(tested_features) > max_features) {
    message2("Subsetting tested features for smoke test: %d -> %d",
             length(tested_features), max_features)
    tested_features <- tested_features[seq_len(max_features)]
  }

  message2("Testing %d features for environment convergence [run_tag=%s]",
           length(tested_features), run_tag)
  pheno <- load_pheno(cohort, tissue)
  cell_props <- load_cell_props(tissue)
  rse <- load_rse(tissue, modality)
  rse <- rse[intersect(tested_features, rownames(rse)), , drop = FALSE]
  mat <- get_assay_matrix(rse, modality)
  colnames(mat) <- sample_ids_from_rse(rse)
  meta <- metadata_for_residualization(rse, pheno, cell_props)
  sample_ids <- intersect(colnames(mat), meta$sample_id)
  mat <- mat[, sample_ids, drop = FALSE]
  meta <- meta[match(sample_ids, meta$sample_id), , drop = FALSE]

  env_pheno <- pheno |>
    filter(brnum %in% sample_ids) |>
    mutate(across(where(is.character), as.factor))
  env_vars <- intersect(get_discovery_env_vars(), colnames(env_pheno))
  if (length(env_vars) == 0) stop("No discovery environmental variables found.")

  fit_env_features <- function(env) {
    samples <- Reduce(intersect, list(colnames(mat), meta$sample_id, env_pheno$brnum))
    dat <- meta[match(samples, meta$sample_id), , drop = FALSE]
    env_dat <- env_pheno[match(samples, env_pheno$brnum), , drop = FALSE]
    dat$env_value <- env_dat[[env]]
    keep <- complete.cases(dat) & !is.na(dat$env_value)
    dat <- dat[keep, , drop = FALSE]
    if (nrow(dat) < 20 || length(unique(dat$env_value)) < 2) {
      return(tibble())
    }
    if (is.character(dat$env_value)) dat$env_value <- as.factor(dat$env_value)
    rhs <- c(
      "env_value",
      setdiff(colnames(dat), c("sample_id", "env_value"))
    )
    form <- as.formula(paste(
      "~",
      paste(vapply(rhs, quote_formula_sym, character(1)), collapse = " + ")
    ))
    design <- model.matrix(form, data = dat)
    if (qr(design)$rank < ncol(design)) return(tibble())

    env_cols <- which(attr(design, "assign") == 1L)
    if (length(env_cols) == 0) return(tibble())

    y <- mat[, dat$sample_id, drop = FALSE]
    y[!is.finite(y)] <- NA_real_
    n_by_feature <- rowSums(!is.na(y))

    fit <- tryCatch(limma::lmFit(y, design), error = function(e) NULL)
    if (is.null(fit)) return(tibble())

    estimate <- fit$coefficients[, env_cols, drop = FALSE]
    std_error <- fit$stdev.unscaled[, env_cols, drop = FALSE] * fit$sigma
    statistic <- estimate / std_error
    p_value <- matrix(
      NA_real_,
      nrow = nrow(statistic),
      ncol = ncol(statistic),
      dimnames = dimnames(statistic)
    )
    for (j in seq_len(ncol(statistic))) {
      p_value[, j] <- 2 * pt(
        abs(statistic[, j]),
        df = fit$df.residual,
        lower.tail = FALSE
      )
    }

    rbindlist(lapply(seq_along(env_cols), function(j) {
      tibble(
        term = colnames(design)[env_cols[[j]]],
        estimate = estimate[, j],
        std.error = std_error[, j],
        statistic = statistic[, j],
        p.value = p_value[, j],
        feature_id = rownames(y),
        env = env,
        n = n_by_feature
      )
    })) |>
      filter(n >= 20, is.finite(p.value))
  }

  env_results <- list()
  idx <- 1L
  for (env in env_vars) {
    message2("Testing environment variable: %s [run_tag=%s]", env, run_tag)
    env_df <- fit_env_features(env)
    if (nrow(env_df) > 0) {
      env_df <- env_df |>
        mutate(
          fdr = p.adjust(p.value, method = "fdr"),
          sig_fdr_10 = !is.na(fdr) & fdr < 0.10,
          sig_p_05 = !is.na(p.value) & p.value < 0.05
        )
      env_results[[idx]] <- env_df
      idx <- idx + 1L
    }
  }

  env_assoc <- bind_rows(env_results) |>
    mutate(
      cohort = cohort, tissue = tissue, population = population,
      vmr_set = vmr_set, modality = modality, run_tag = run_tag
    )
  safe_fwrite(env_assoc, file.path(out_dir, "feature_environment_associations.tsv.gz"),
              sep = "\t")

  ## Build feature context table: classify each feature as linked exclusively to
  ## heritable VMRs, non-heritable VMRs, or both, based on sig. DNAm associations
  ## (FDR < 0.10 in the local methylation-feature model).
  dnam_linked <- assoc |>
    filter(sig_fdr_10, h2_category %in% c("Heritable", "Non-heritable")) |>
    mutate(group = as.character(h2_category)) |>
    distinct(feature_id, group)

  feature_group <- dnam_linked |>
    group_by(feature_id) |>
    summarise(
      linked_heritable = any(group == "Heritable"),
      linked_nonheritable = any(group == "Non-heritable"),
      feature_context = case_when(
        linked_nonheritable & !linked_heritable ~ "Non-heritable only",
        linked_heritable & !linked_nonheritable ~ "Heritable only",
        linked_nonheritable & linked_heritable ~ "Both",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    ) |>
    filter(feature_context %in% c("Non-heritable only", "Heritable only"))

  run_fisher <- function(env, term) {
    env_sig <- env_assoc |>
      filter(env == !!env, term == !!term) |>
      transmute(feature_id, env_linked = sig_p_05)
    tab_df <- feature_group |>
      inner_join(env_sig, by = "feature_id")
    if (nrow(tab_df) == 0 || length(unique(tab_df$env_linked)) < 2) {
      return(tibble())
    }
    tab <- table(tab_df$env_linked, tab_df$feature_context)
    needed <- c("Non-heritable only", "Heritable only")
    if (!all(needed %in% colnames(tab)) || nrow(tab) < 2) return(tibble())
    ft <- fisher.test(tab[, needed])
    tibble(
      env = env,
      term = term,
      odds_ratio = unname(ft$estimate),
      ci_lo = ft$conf.int[[1]],
      ci_hi = ft$conf.int[[2]],
      p_value = ft$p.value,
      n_nonheritable_features = sum(tab_df$feature_context == "Non-heritable only"),
      n_heritable_features = sum(tab_df$feature_context == "Heritable only"),
      n_env_linked_nonheritable = sum(tab_df$env_linked &
                                        tab_df$feature_context == "Non-heritable only"),
      n_env_linked_heritable = sum(tab_df$env_linked &
                                     tab_df$feature_context == "Heritable only")
    )
  }

  pairs <- env_assoc |> distinct(env, term)
  enrichment <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
    run_fisher(pairs$env[[i]], pairs$term[[i]])
  }))
  if (nrow(enrichment) > 0) {
    enrichment <- enrichment |>
      mutate(
        fdr = p.adjust(p_value, method = "fdr"),
        cohort = cohort, tissue = tissue, population = population,
        vmr_set = vmr_set, modality = modality, run_tag = run_tag
      )
  }
  safe_fwrite(enrichment, file.path(out_dir, "environment_convergence_enrichment.tsv"),
              sep = "\t")

  ## Rank-based Wilcoxon test (no threshold)
  run_wilcoxon_rank <- function(env, term) {
    env_scores <- env_assoc |>
      filter(env == !!env, term == !!term, is.finite(p.value)) |>
      transmute(feature_id, neg_log10_p = -log10(p.value + .Machine$double.eps))
    tab_df <- feature_group |>
      inner_join(env_scores, by = "feature_id")
    if (nrow(tab_df) == 0) return(tibble())
    contexts <- unique(tab_df$feature_context)
    if (!all(c("Heritable only", "Non-heritable only") %in% contexts)) return(tibble())
    x_h  <- tab_df$neg_log10_p[tab_df$feature_context == "Heritable only"]
    x_nh <- tab_df$neg_log10_p[tab_df$feature_context == "Non-heritable only"]
    if (length(x_h) < 3 || length(x_nh) < 3) return(tibble())
    wt <- wilcox.test(x_h, x_nh, exact = FALSE, conf.int = TRUE)
    tibble(
      env = env,
      term = term,
      statistic = unname(wt$statistic),
      p_value = wt$p.value,
      location_shift = unname(wt$estimate),
      ci_lo = wt$conf.int[[1]],
      ci_hi = wt$conf.int[[2]],
      median_neg_log10p_heritable = median(x_h),
      median_neg_log10p_nonheritable = median(x_nh),
      n_heritable_features = length(x_h),
      n_nonheritable_features = length(x_nh)
    )
  }

  wilcoxon_enrichment <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
    run_wilcoxon_rank(pairs$env[[i]], pairs$term[[i]])
  }))
  if (nrow(wilcoxon_enrichment) > 0) {
    wilcoxon_enrichment <- wilcoxon_enrichment |>
      mutate(
        fdr = p.adjust(p_value, method = "fdr"),
        cohort = cohort, tissue = tissue, population = population,
        vmr_set = vmr_set, modality = modality, run_tag = run_tag
      )
  }
  safe_fwrite(wilcoxon_enrichment,
              file.path(out_dir, "environment_convergence_wilcoxon.tsv"),
              sep = "\t")

  ## VMR-level aggregation via Cauchy combination
  cauchy_combine <- function(p_vals) {
    p_vals <- p_vals[is.finite(p_vals) & p_vals > 0 & p_vals <= 1]
    if (length(p_vals) == 0) return(NA_real_)
    t_stat <- mean(tan((0.5 - p_vals) * pi))
    pt(t_stat, df = 1, lower.tail = FALSE)
  }

  run_vmr_level <- function(env, term) {
    env_p <- env_assoc |>
      filter(env == !!env, term == !!term, is.finite(p.value)) |>
      select(feature_id, p.value)

    vmr_scores <- assoc |>
      filter(h2_category %in% c("Heritable", "Non-heritable"), sig_fdr_10) |>
      select(vmr_id, h2_category, feature_id) |>
      distinct() |>
      inner_join(env_p, by = "feature_id") |>
      group_by(vmr_id, h2_category) |>
      summarise(
        cauchy_p = cauchy_combine(p.value),
        n_features = n(),
        .groups = "drop"
      ) |>
      filter(is.finite(cauchy_p))

    if (nrow(vmr_scores) == 0) return(tibble())
    contexts <- unique(vmr_scores$h2_category)
    if (!all(c("Heritable", "Non-heritable") %in% contexts)) return(tibble())

    x_h  <- -log10(vmr_scores$cauchy_p[vmr_scores$h2_category == "Heritable"] +
                     .Machine$double.eps)
    x_nh <- -log10(vmr_scores$cauchy_p[vmr_scores$h2_category == "Non-heritable"] +
                     .Machine$double.eps)
    if (length(x_h) < 3 || length(x_nh) < 3) return(tibble())

    wt <- wilcox.test(x_h, x_nh, exact = FALSE, conf.int = TRUE)
    tibble(
      env = env,
      term = term,
      statistic = unname(wt$statistic),
      p_value = wt$p.value,
      location_shift = unname(wt$estimate),
      ci_lo = wt$conf.int[[1]],
      ci_hi = wt$conf.int[[2]],
      median_vmr_score_heritable = median(x_h),
      median_vmr_score_nonheritable = median(x_nh),
      n_heritable_vmrs = length(x_h),
      n_nonheritable_vmrs = length(x_nh)
    )
  }

  vmr_level_enrichment <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
    run_vmr_level(pairs$env[[i]], pairs$term[[i]])
  }))
  if (nrow(vmr_level_enrichment) > 0) {
    vmr_level_enrichment <- vmr_level_enrichment |>
      mutate(
        fdr = p.adjust(p_value, method = "fdr"),
        cohort = cohort, tissue = tissue, population = population,
        vmr_set = vmr_set, modality = modality, run_tag = run_tag
      )
  }
  safe_fwrite(vmr_level_enrichment,
              file.path(out_dir, "environment_convergence_vmr_level.tsv"),
              sep = "\t")

  message2("Saved environment convergence outputs to %s [run_tag=%s]", out_dir, run_tag)
}

#### Reproducibility ####
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
