#### Nearest gene/PSI proximity by VMR genetic architecture ####
#
# Genome-space distances from each elastic-net VMR to the nearest annotated gene
# body and nearest PSI feature (inputs/counts). Independent of ABC vs nearest-gene
# association outputs under regulatory_context/_m/.../expression/{abc,nearest_gene_*}/.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(GenomicRanges)
  library(broom)
})

source(here("heritability", "elastic_net_model", "all_individuals",
            "tissue_comparison", "regulatory_context", "_h",
            "00.regulatory_context_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
cohort <- ifelse(length(args) >= 1, args[[1]], "all_individuals")
tissue <- ifelse(length(args) >= 2, tolower(args[[2]]), "dlpfc")
population <- ifelse(length(args) >= 3, toupper(args[[3]]), "AA")
vmr_set <- ifelse(length(args) >= 4, args[[4]], "shared")
vmr_set <- validate_vmr_set(cohort, population, vmr_set)

if (should_skip_shared_duplicate_population(population, vmr_set)) {
  message2(
    paste0(
      "Skipping duplicate shared-VMR proximity run for population=%s ",
      "(identical genomic distances and h2_category; see %s)."
    ),
    population,
    regctx_output_dir(cohort, tissue, SHARED_VMR_CANONICAL_POPULATION,
                      "proximity", vmr_set = vmr_set)
  )
  quit(save = "no", status = 0)
}

out_dir <- regctx_output_dir(cohort, tissue, population, "proximity",
                             vmr_set = vmr_set)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nearest_feature <- function(vmr_df, annot_df, feature_type) {
  vmr_gr <- GRanges(vmr_df$seqnames, IRanges(vmr_df$start, vmr_df$end))
  feat_gr <- GRanges(annot_df$chrom, IRanges(annot_df$start, annot_df$end))
  idx <- nearest(vmr_gr, feat_gr, ignore.strand = TRUE)
  tibble(
    vmr_id = vmr_df$vmr_id,
    feature_type = feature_type,
    nearest_feature_id = annot_df$feature_id[idx],
    nearest_feature_name = annot_df$feature_name[idx],
    nearest_distance = distance(vmr_gr, feat_gr[idx])
  )
}

enet <- load_enet(cohort, tissue, population, vmr_set) |>
  mutate(
    vmr_id = paste(seqnames, start, end, sep = "_"),
    vmr_length = end - start,
    h2_category = factor(h2_category, levels = H2_GROUP_LEVELS)
  )

gene_annot <- safe_read(here("inputs", "counts", "gene-annotation.tsv")) |>
  transmute(
    chrom = normalize_chr(chrom),
    start = as.integer(start),
    end = as.integer(end),
    feature_id = strip_ensembl_version(gene_id),
    feature_name = gene_name
  ) |>
  filter(!is.na(chrom), !is.na(start), !is.na(end), end >= start)

psi_annot <- safe_read(here("inputs", "counts", "psi-annotation.tsv")) |>
  transmute(
    chrom = normalize_chr(chrom),
    start = as.integer(start),
    end = as.integer(end),
    feature_id = psi_uid,
    feature_name = paste0(gene_name, "|", event_type)
  ) |>
  filter(!is.na(chrom), !is.na(start), !is.na(end), end >= start)

proximity <- bind_rows(
  nearest_feature(enet, gene_annot, "gene"),
  nearest_feature(enet, psi_annot, "psi")
) |>
  left_join(enet, by = "vmr_id") |>
  mutate(
    log_nearest_distance = log1p(nearest_distance),
    within_10kb = nearest_distance <= 10000,
    within_100kb = nearest_distance <= 100000
  )

safe_fwrite(proximity, file.path(out_dir, "vmr_nearest_feature_distances.tsv"),
            sep = "\t")

make_bins <- function(x, n_bins, prefix) {
  br <- unique(quantile(x, probs = seq(0, 1, length.out = n_bins + 1),
                        na.rm = TRUE))
  if (length(br) < 2) return(factor(rep(NA_character_, length(x))))
  cut(x, breaks = br, include.lowest = TRUE,
      labels = paste0(prefix, seq_len(length(br) - 1)))
}

prox_high_conf <- proximity |>
  filter(r_squared_cv > 0.3, !is.na(h2_unscaled)) |>
  group_by(feature_type) |>
  mutate(
    h2_tertile = make_bins(h2_unscaled, 3, "T"),
    h2_quintile = make_bins(h2_unscaled, 5, "Q")
  ) |>
  ungroup()

quintile_summary <- prox_high_conf |>
  pivot_longer(c(h2_tertile, h2_quintile),
               names_to = "bin_type", values_to = "h2_bin") |>
  filter(!is.na(h2_bin)) |>
  group_by(cohort = cohort, tissue = tissue, population = population,
           vmr_set = vmr_set,
           feature_type, bin_type, h2_bin) |>
  summarise(
    n_vmrs = n(),
    median_distance = median(nearest_distance, na.rm = TRUE),
    q25_distance = quantile(nearest_distance, 0.25, na.rm = TRUE),
    q75_distance = quantile(nearest_distance, 0.75, na.rm = TRUE),
    prop_within_10kb = mean(within_10kb, na.rm = TRUE),
    prop_within_100kb = mean(within_100kb, na.rm = TRUE),
    .groups = "drop"
  )
safe_fwrite(quintile_summary, file.path(out_dir, "proximity_h2_bin_summary.tsv"),
            sep = "\t")

group_summary <- proximity |>
  filter(h2_category %in% c("Heritable", "Non-heritable", "Low prediction")) |>
  group_by(cohort = cohort, tissue = tissue, population = population,
           vmr_set = vmr_set,
           feature_type, h2_category) |>
  summarise(
    n_vmrs = n(),
    median_distance = median(nearest_distance, na.rm = TRUE),
    q25_distance = quantile(nearest_distance, 0.25, na.rm = TRUE),
    q75_distance = quantile(nearest_distance, 0.75, na.rm = TRUE),
    prop_within_10kb = mean(within_10kb, na.rm = TRUE),
    prop_within_100kb = mean(within_100kb, na.rm = TRUE),
    .groups = "drop"
  )
safe_fwrite(group_summary, file.path(out_dir, "proximity_h2_group_summary.tsv"),
            sep = "\t")

test_rows <- list()
for (ft in c("gene", "psi")) {
  df_ft <- proximity |>
    filter(feature_type == ft,
           h2_category %in% c("Heritable", "Non-heritable"),
           !is.na(log_nearest_distance))

  if (nrow(df_ft) > 10) {
    wt <- wilcox.test(log_nearest_distance ~ h2_category, data = df_ft)
    test_rows[[length(test_rows) + 1]] <- tibble(
      cohort = cohort, tissue = tissue, population = population,
      vmr_set = vmr_set,
      feature_type = ft, comparison = "Heritable_vs_non-heritable",
      outcome = "log_nearest_distance",
      test = "wilcoxon",
      statistic = unname(wt$statistic),
      p_value = wt$p.value
    )

    for (threshold in c("within_10kb", "within_100kb")) {
      if (length(unique(df_ft[[threshold]])) > 1) {
        fit <- glm(
          as.formula(paste0(threshold,
                            " ~ h2_category + log(vmr_length) + log(num_snps + 1)")),
          data = df_ft,
          family = binomial
        )
        test_rows[[length(test_rows) + 1]] <- tidy(fit, conf.int = TRUE,
                                                   exponentiate = TRUE) |>
          filter(grepl("^h2_category", term)) |>
          transmute(
            cohort = cohort, tissue = tissue, population = population,
            vmr_set = vmr_set,
            feature_type = ft,
            comparison = "Heritable_vs_non-heritable_adjusted",
            outcome = threshold,
            test = "logistic",
            term = term,
            estimate = estimate,
            conf.low = conf.low,
            conf.high = conf.high,
            p_value = p.value
          )
      }
    }
  }

  df_q <- prox_high_conf |>
    filter(feature_type == ft, !is.na(h2_quintile), !is.na(log_nearest_distance))
  if (nrow(df_q) > 10 && length(unique(df_q$h2_quintile)) > 1) {
    kt <- kruskal.test(log_nearest_distance ~ h2_quintile, data = df_q)
    test_rows[[length(test_rows) + 1]] <- tibble(
      cohort = cohort, tissue = tissue, population = population,
      vmr_set = vmr_set,
      feature_type = ft, comparison = "h2_quintile",
      outcome = "log_nearest_distance",
      test = "kruskal",
      statistic = unname(kt$statistic),
      p_value = kt$p.value
    )
  }
}

tests <- bind_rows(test_rows)
if (nrow(tests) > 0) {
  tests <- tests |> mutate(fdr = p.adjust(p_value, method = "fdr"))
}
safe_fwrite(tests, file.path(out_dir, "proximity_group_tests.tsv"), sep = "\t")

message2("Saved proximity outputs to %s", out_dir)

#### Reproducibility ####
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
