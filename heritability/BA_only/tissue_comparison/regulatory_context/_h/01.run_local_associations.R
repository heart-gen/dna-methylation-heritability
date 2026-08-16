#### Local VMR methylation - expression (ABC + nearest gene) / PSI associations ####

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

## -----------------------------------------------------------------------------
## One association pass (shared by ABC, nearest gene, and PSI)
## -----------------------------------------------------------------------------

run_association_pass <- function(
    links,
    mat_full,
    meta,
    sample_ids,
    pheno,
    cohort,
    tissue,
    population,
    modality_label,
    assay_name,
    meth,
    out_dir,
    max_pairs,
    n_rna_samples_total,
    n_features_total
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  varying_ids <- feature_ids_with_variation(mat_full, sample_ids)
  n_pairs_before_constant_filter <- nrow(links)
  links <- links |> filter(feature_id %in% varying_ids)
  n_pairs_after_constant_filter <- nrow(links)
  message2(
    "%s: assay rows with variation %d / %d; VMR-feature links %d -> %d after dropping constant / all-zero features",
    modality_label, length(varying_ids), nrow(mat_full),
    n_pairs_before_constant_filter, n_pairs_after_constant_filter
  )
  if (nrow(links) == 0) {
    stop("No VMR-feature links remain after removing constant or all-zero features")
  }

  if (!is.na(max_pairs) && nrow(links) > max_pairs) {
    message2("Subsetting links for smoke test: %d -> %d", nrow(links), max_pairs)
    links <- links |> slice_head(n = max_pairs)
  }

  mat <- mat_full[intersect(unique(links$feature_id), rownames(mat_full)),
                  sample_ids, drop = FALSE]
  meta <- meta[match(sample_ids, meta$sample_id), , drop = FALSE]

  coverage <- tibble(
    cohort = cohort, tissue = tissue, modality = modality_label, population = population,
    assay = assay_name, n_control_pheno = nrow(pheno),
    n_rna_samples = n_rna_samples_total,
    n_model_samples_pre_pair_filter = length(sample_ids), n_features_total = n_features_total,
    n_features_assay_varying = length(varying_ids), n_pairs_before_constant_filter = n_pairs_before_constant_filter,
    n_pairs_after_constant_filter = n_pairs_after_constant_filter, n_features_linked = nrow(mat),
    model = "feature ~ VMR_methylation + RNA_covariates"
  )
  safe_fwrite(coverage, file.path(out_dir, "sample_feature_coverage.tsv"), sep = "\t")

  safe_fwrite(links, file.path(out_dir, "tested_vmr_feature_links.tsv.gz"),
              sep = "\t")

  tested_vmrs <- unique(links$vmr_id)
  meth_wide <- meth |>
    filter(vmr_id %in% tested_vmrs) |>
    dplyr::select(brnum, vmr_id, meth) |>
    distinct() |>
    tidyr::pivot_wider(names_from = "vmr_id", values_from = "meth") |>
    as.data.frame()
  rownames(meth_wide) <- meth_wide$brnum
  meth_wide$brnum <- NULL

  common_samples <- Reduce(intersect, list(rownames(meth_wide), colnames(mat),
                                          meta$sample_id))
  meta_common <- meta[match(common_samples, meta$sample_id), , drop = FALSE]
  cov_cols <- setdiff(colnames(meta_common), "sample_id")

  form_x <- as.formula(paste(
    "~",
    paste(vapply(c("meth", cov_cols), quote_formula_sym, character(1)),
          collapse = " + ")
  ))

  na_fit <- function(n = 0L) {
    tibble(beta = NA_real_, se = NA_real_, t = NA_real_,
           p = NA_real_, n = n, r = NA_real_)
  }

  fit_pair_single <- function(vmr_id, feature_id) {
    if (!vmr_id %in% colnames(meth_wide) || !feature_id %in% rownames(mat)) {
      return(na_fit(0L))
    }
    dat      <- meta_common
    dat$meth <- as.numeric(meth_wide[common_samples, vmr_id])
    dat$feature <- as.numeric(mat[feature_id, common_samples])
    covars   <- setdiff(colnames(dat), "sample_id")
    keep     <- complete.cases(dat[, covars, drop = FALSE]) &
      is.finite(dat$meth) & is.finite(dat$feature)
    dat <- dat[keep, , drop = FALSE]
    n <- nrow(dat)
    if (n < 20 || length(unique(dat$meth)) < 3) {
      return(na_fit(n))
    }
    rhs <- c("meth", setdiff(colnames(dat), c("sample_id", "feature", "meth")))
    form <- as.formula(paste(
      "feature ~",
      paste(vapply(rhs, quote_formula_sym, character(1)), collapse = " + ")
    ))
    design <- model.matrix(form, data = dat)
    if (qr(design)$rank < ncol(design)) {
      return(na_fit(n))
    }
    fit <- lm(form, data = dat)
    sm  <- summary(fit)$coefficients
    tibble(
      beta = sm["meth", "Estimate"], se = sm["meth", "Std. Error"],
      t = sm["meth", "t value"], p = sm["meth", "Pr(>|t|)"],
      n = n, r = suppressWarnings(cor(dat$meth, dat$feature))
    )
  }

  infer_meth_xtinv <- function(X, XtX_inv, y, meth_j, n, p, meth_v, feat_v) {
    if (anyNA(y)) {
      return(na_fit(as.integer(n)))
    }
    df_r <- n - p
    if (df_r <= 0L) {
      return(na_fit(as.integer(n)))
    }
    coef   <- as.vector(XtX_inv %*% crossprod(X, y))
    fit    <- as.vector(X %*% coef)
    res    <- y - fit
    sigma2 <- sum(res^2) / df_r
    v_m    <- XtX_inv[meth_j, meth_j] * sigma2
    if (!is.finite(v_m) || v_m <= 0) {
      return(na_fit(as.integer(n)))
    }
    beta_m <- coef[meth_j]
    se_m   <- sqrt(v_m)
    if (!is.finite(beta_m) || is.na(beta_m)) {
      return(na_fit(as.integer(n)))
    }
    tstat <- beta_m / se_m
    pval <- 2 * pt(-abs(tstat), df = df_r)
    tibble(
      beta = beta_m, se = se_m, t = tstat, p = pval, n = as.integer(n),
      r = suppressWarnings(cor(meth_v, feat_v))
    )
  }

  message2("Testing %d VMR-feature pairs (batched by VMR) [%s]",
           nrow(links), modality_label)
  results   <- vector("list", nrow(links))
  pair_done <- 0L
  vmr_ids   <- unique(links$vmr_id)
  n_vmrs    <- length(vmr_ids)

  for (v in seq_len(n_vmrs)) {
    vmr_id <- vmr_ids[[v]]
    idx    <- which(links$vmr_id == vmr_id)
    if (length(idx) == 0L) next

    if (!vmr_id %in% colnames(meth_wide)) {
      for (j in idx) results[[j]] <- na_fit(0L)
      pair_done <- pair_done + length(idx)
      next
    }

    meth_vec  <- as.numeric(meth_wide[common_samples, vmr_id])
    dat0      <- meta_common
    dat0$meth <- meth_vec
    base_keep <- stats::complete.cases(dat0[, cov_cols, drop = FALSE]) &
      is.finite(dat0$meth)
    n_base <- sum(base_keep)
    meth_b <- meth_vec[base_keep]

    if (n_base < 20L || length(unique(meth_b)) < 3L) {
      for (j in idx) results[[j]] <- na_fit(as.integer(n_base))
      pair_done <- pair_done + length(idx)
      next
    }

    dat_base <- dat0[base_keep, , drop = FALSE]
    X      <- stats::model.matrix(form_x, data = dat_base)
    p_col  <- ncol(X)
    n_row  <- nrow(X)
    meth_j <- match("meth", colnames(X))
    if (is.na(meth_j)) {
      for (j in idx) results[[j]] <- fit_pair_single(vmr_id, links$feature_id[[j]])
      pair_done <- pair_done + length(idx)
      next
    }

    XtX       <- crossprod(X)
    XtX_inv   <- tryCatch(solve(XtX), error = function(e) NULL)
    use_batch <- !is.null(XtX_inv) && (n_row > p_col)

    if (use_batch) {
      for (j in idx) {
        fid <- links$feature_id[[j]]
        if (!fid %in% rownames(mat)) {
          results[[j]] <- na_fit(0L)
          next
        }
        y_full <- as.numeric(mat[fid, common_samples])
        if (!all(is.finite(y_full[base_keep]))) {
          results[[j]] <- fit_pair_single(vmr_id, fid)
          next
        }
        y <- y_full[base_keep]
        if (stats::var(y) <= .Machine$double.eps) {
          results[[j]] <- na_fit(as.integer(n_base))
          next
        }
        results[[j]] <- infer_meth_xtinv(X, XtX_inv, y, meth_j, n_row, p_col,
                                         meth_b, y)
      }
    } else {
      for (j in idx) {
        results[[j]] <- fit_pair_single(vmr_id, links$feature_id[[j]])
      }
    }

    pair_done <- pair_done + length(idx)
    if (pair_done %% 50000L == 0L || v %% max(1L, floor(n_vmrs / 20L)) == 0L) {
      message2("  tested %d / %d pairs (%d / %d VMRs) [%s]",
               pair_done, nrow(links), v, n_vmrs, modality_label)
    }
  }

  assoc <- bind_cols(links, bind_rows(results)) |>
    mutate(
      fdr = p.adjust(p, method = "fdr"),
      abs_beta = abs(beta),
      sig_fdr_05 = !is.na(fdr) & fdr < 0.05,
      sig_fdr_10 = !is.na(fdr) & fdr < 0.10
    )

  safe_fwrite(assoc, file.path(out_dir, "vmr_feature_associations.tsv.gz"),
              sep = "\t")

  vmr_summary <- assoc |>
    group_by(cohort = cohort, tissue, population, modality = modality_label, vmr_id,
             seqnames, start, end, h2_category, h2_unscaled) |>
    summarise(
      n_pairs_tested = sum(!is.na(p)),
      n_features_tested = n_distinct(feature_id[!is.na(p)]),
      n_sig_fdr_05 = sum(sig_fdr_05, na.rm = TRUE),
      n_sig_fdr_10 = sum(sig_fdr_10, na.rm = TRUE),
      min_p = suppressWarnings(min(p, na.rm = TRUE)),
      min_fdr = suppressWarnings(min(fdr, na.rm = TRUE)),
      max_abs_beta = suppressWarnings(max(abs_beta, na.rm = TRUE)),
      strongest_beta = beta[which.max(replace(abs_beta, is.na(abs_beta), -Inf))][1],
      .groups = "drop"
    ) |>
    mutate(
      min_p = ifelse(is.infinite(min_p), NA_real_, min_p),
      min_fdr = ifelse(is.infinite(min_fdr), NA_real_, min_fdr),
      max_abs_beta = ifelse(is.infinite(max_abs_beta), NA_real_, max_abs_beta)
    )

  safe_fwrite(vmr_summary, file.path(out_dir, "vmr_association_summary.tsv"),
              sep = "\t")

  message2("Saved [%s] outputs to %s", modality_label, out_dir)
}

args       <- commandArgs(trailingOnly = TRUE)
cohort     <- ifelse(length(args) >= 1, args[[1]], "BA_only")
tissue     <- ifelse(length(args) >= 2, tolower(args[[2]]), "dlpfc")
modality   <- ifelse(length(args) >= 3, tolower(args[[3]]), "expression")
population <- ifelse(length(args) >= 4, toupper(args[[4]]), "AA")
window     <- ifelse(length(args) >= 5, as.integer(args[[5]]), 250000L)
max_pairs  <- ifelse(length(args) >= 6, as.integer(args[[6]]), NA_integer_)

if (!cohort %in% c("BA_only", "all_individuals")) {
  stop("cohort must be BA_only or all_individuals")
}
if (!tissue %in% TISSUES) stop("Unknown tissue: ", tissue)
if (!modality %in% c("expression", "psi")) {
  stop("modality must be expression or psi")
}
if (!population %in% c("AA", "EA")) stop("population must be AA or EA")
if (cohort == "BA_only" && population != "AA") {
  stop("BA_only supports population AA only")
}

message2("Running local associations: cohort=%s tissue=%s modality=%s population=%s",
         cohort, tissue, modality, population)

## Load RNA/PSI features and covariates.
pheno      <- load_pheno(cohort, tissue)
cell_props <- load_cell_props(tissue)
rse        <- load_rse(tissue, modality)
mat_full   <- get_assay_matrix(rse, modality)
assay_name <- attr(mat_full, "assay_name")
n_rna_samples_total <- ncol(mat_full)
n_features_total    <- nrow(mat_full)
colnames(mat_full) <- sample_ids_from_rse(rse)
feature_annot <- feature_annotation(rse, modality)

meta       <- metadata_for_residualization(rse, pheno, cell_props)
sample_ids <- intersect(colnames(mat_full), meta$sample_id)

if (modality == "expression") {
  out_dir_abc <- here(
    "heritability", "elastic_net_model", cohort,
    "tissue_comparison", "regulatory_context", "_m",
    tissue, population, "expression", "abc"
  )
  links_abc <- load_abc_links(cohort, population) |>
    filter(tissue == !!tissue) |>
    inner_join(feature_annot,
               by = c("TargetGene_base" = "gene_id_base")) |>
    transmute(
      seqnames, start, end, vmr_id, h2_category, h2_unscaled,
      tissue, population = population,
      feature_id, feature_label = coalesce(gene_name, TargetGene_name),
      target_id = TargetGene, target_name = TargetGene_name,
      link_type = "ABC",
      link_score = ABC.Score,
      distance = distance
    ) |>
    distinct()

  links_ng <- tryCatch(
    load_nearest_gene_links(cohort, tissue, population, window) |>
      inner_join(feature_annot, by = "gene_id_base") |>
      transmute(
        seqnames, start, end, vmr_id, h2_category, h2_unscaled,
        tissue, population, feature_id,
        feature_label = coalesce(gene_name, nearest_gene_symbol_within_250kb),
        target_id = gene_id,
        target_name = coalesce(gene_name, nearest_gene_symbol_within_250kb),
        link_type = paste0("nearest_gene_", window / 1000, "kb"),
        link_score = NA_real_,
        distance = distance_to_nearest_gene_within_250kb
      ) |>
      distinct(),
    error = function(e) {
      message2("Nearest-gene links unavailable (%s); ABC-only run.",
               conditionMessage(e))
      tibble()
    }
  )

  if (nrow(links_abc) == 0) stop("No ABC VMR–gene links for this run.")

  mt   <- load_methylation(cohort, tissue, pheno)
  meth <- mt |> filter(brnum %in% sample_ids)
  if (nrow(meth) == 0) stop("No methylation rows overlap RNA/covariate samples.")

  run_association_pass(
    links = links_abc,
    mat_full = mat_full,
    meta = meta,
    sample_ids = sample_ids,
    pheno = pheno,
    cohort = cohort,
    tissue = tissue,
    population = population,
    modality_label = "expression",
    assay_name = assay_name,
    meth = meth,
    out_dir = out_dir_abc,
    max_pairs = max_pairs,
    n_rna_samples_total = n_rna_samples_total,
    n_features_total = n_features_total
  )

  if (nrow(links_ng) > 0) {
    out_dir_ng <- here(
      "heritability", "elastic_net_model", cohort,
      "tissue_comparison", "regulatory_context", "_m",
      tissue, population, "expression",
      paste0("nearest_gene_window_", window / 1000, "kb")
    )
    message2("Running nearest-gene association pass (%d links)", nrow(links_ng))
    run_association_pass(
      links = links_ng,
      mat_full = mat_full,
      meta = meta,
      sample_ids = sample_ids,
      pheno = pheno,
      cohort = cohort,
      tissue = tissue,
      population = population,
      modality_label = "expression_nearest_gene",
      assay_name = assay_name,
      meth = meth,
      out_dir = out_dir_ng,
      max_pairs = max_pairs,
      n_rna_samples_total = n_rna_samples_total,
      n_features_total = n_features_total
    )
  } else {
    message2("Skipping nearest-gene association pass (no links after filters).")
  }

  message2("Expression passes finished (ABC + nearest gene when available).")
} else {
  out_dir <- here("heritability", "elastic_net_model", cohort,
                  "tissue_comparison", "regulatory_context", "_m",
                  tissue, population, modality,
                  paste0("window_", window / 1000, "kb"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  links <- load_psi_links(cohort, tissue, population, window) |>
    inner_join(feature_annot, by = c("psi_uid" = "psi_uid")) |>
    transmute(
      seqnames, start, end, vmr_id, h2_category, h2_unscaled,
      tissue, population, feature_id,
      feature_label = paste0(psi_uid, "|", gene_name, "|", event_type),
      target_id = gene_id, target_name = gene_name,
      link_type = paste0("local_", window / 1000, "kb"),
      link_score = NA_real_,
      distance = distance
    ) |>
    distinct()

  if (nrow(links) == 0) stop("No local feature links available for this run.")

  mt   <- load_methylation(cohort, tissue, pheno)
  meth <- mt |> filter(brnum %in% sample_ids)
  if (nrow(meth) == 0) stop("No methylation rows overlap RNA/covariate samples.")

  run_association_pass(
    links = links,
    mat_full = mat_full,
    meta = meta,
    sample_ids = sample_ids,
    pheno = pheno,
    cohort = cohort,
    tissue = tissue,
    population = population,
    modality_label = "psi",
    assay_name = assay_name,
    meth = meth,
    out_dir = out_dir,
    max_pairs = max_pairs,
    n_rna_samples_total = n_rna_samples_total,
    n_features_total = n_features_total
  )

  message2("Saved outputs to %s", out_dir)
}

#### Reproducibility ####
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
