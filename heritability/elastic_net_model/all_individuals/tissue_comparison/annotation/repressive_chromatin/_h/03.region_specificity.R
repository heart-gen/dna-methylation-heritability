#### Region-Specificity and Sharing of Repressive Heritable VMRs ####
##
## Extends the repressive chromatin enrichment analysis (01.) with a
## cross-tissue sharing dimension, asking:
##
##   Q1. Is the repressive heritable compartment shared or region-specific?
##   Q2. Is repressive enrichment especially strong in caudate-specific VMRs?
##   Q3. Are repressive heritable VMRs more shared than other heritable VMRs?
##   Q4. Are repressive heritable VMRs more shared than non-heritable
##       intergenic VMRs?
##   Q5. Do overlapping repressive heritable VMRs show similar h² across
##       regions?
##
## Uses:
##   - vmr_repressive_overlap.tsv (per-VMR repressive flags, from 01.)
##   - F_0.25 pairwise overlap TSVs (from upset_plot pipeline)
##
## VMRs are independently called per tissue; "sharing" = ≥25% reciprocal
## overlap between tissue-specific VMR calls (F_0.25 threshold).

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

## Configuration

TISSUES <- c("Caudate", "DLPFC", "Hippocampus")
TISSUE_PAIRS <- list(
  c("caudate", "dlpfc"),
  c("caudate", "hippocampus"),
  c("hippocampus", "dlpfc")
)
PAIR_LABELS <- c(
  "caudate_dlpfc"          = "Caudate-DLPFC",
  "caudate_hippocampus"    = "Caudate-Hippocampus",
  "hippocampus_dlpfc"      = "Hippocampus-DLPFC"
)
TISSUE_LABELS <- c(
  caudate = "Caudate",
  dlpfc = "DLPFC",
  hippocampus = "Hippocampus"
)

# Primary repressive definition
REPRESSIVE_COL <- "in_BroadRepressive"
# Sensitivity
SENSITIVITY_COLS <- c("in_H3K9me3", "in_Het", "in_Quies")

## Paths

REPRESSIVE_F <- here("heritability", "elastic_net_model", "BA_only",
                     "tissue_comparison", "annotation", "repressive_chromatin",
                     "_m", "vmr_repressive_overlap.tsv")
OVERLAP_DIR  <- here("heritability", "elastic_net_model", "BA_only",
                     "tissue_comparison", "upset_plot", "_m",
                     "F_0.25", "percent_overlap")
OUT_DIR      <- here("heritability", "elastic_net_model", "BA_only",
                     "tissue_comparison", "annotation", "repressive_chromatin", "_m")

## Functions

ac_ci <- function(n_annot, n, z = qnorm(0.975)) {
  n_tilde <- n + z^2
  p_tilde <- (n_annot + z^2 / 2) / n_tilde
  lo <- pmax(0, p_tilde - z * sqrt(p_tilde * (1 - p_tilde) / n_tilde))
  hi <- pmin(1, p_tilde + z * sqrt(p_tilde * (1 - p_tilde) / n_tilde))
  list(lo = lo, hi = hi)
}

run_fishers_2x2 <- function(n_a_pos, n_a_total, n_b_pos, n_b_total,
                            tissue, question, group_a, group_b) {
  tab <- matrix(
    c(n_a_pos, n_a_total - n_a_pos,
      n_b_pos, n_b_total - n_b_pos),
    nrow = 2, byrow = TRUE,
    dimnames = list(c(group_a, group_b), c("yes", "no"))
  )
  ft <- fisher.test(tab)
  tibble(
    tissue    = tissue,
    question  = question,
    group_a   = group_a,
    group_b   = group_b,
    n_a_pos   = n_a_pos,
    n_a_total = n_a_total,
    rate_a    = n_a_pos / n_a_total,
    n_b_pos   = n_b_pos,
    n_b_total = n_b_total,
    rate_b    = n_b_pos / n_b_total,
    or        = ft$estimate,
    ci_lo     = ft$conf.int[1],
    ci_hi     = ft$conf.int[2],
    p_value   = ft$p.value
  )
}

fisher_z_test <- function(rho1, n1, rho2, n2) {
  z1 <- atanh(rho1)
  z2 <- atanh(rho2)
  se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  z_stat <- (z1 - z2) / se
  p_value <- 2 * pnorm(-abs(z_stat))
  list(z_stat = z_stat, p_value = p_value)
}

## ============================================================
## Step 1: Load repressive overlap data
## ============================================================

cat("Loading repressive overlap data...\n")
repr <- fread(REPRESSIVE_F)
cat(sprintf("  Loaded %d VMRs across %d tissues\n", nrow(repr), length(unique(repr$tissue))))

# Create VMR key for matching
repr[, vmr_key := paste(seqnames, start, end, sep = ":")]

# Define subgroups
repr[, repressive := as.integer(get(REPRESSIVE_COL) == 1)]
repr[, subgroup := fcase(
  h2_category == "Heritable" & repressive == 1, "repressive_heritable",
  h2_category == "Heritable" & repressive == 0, "nonrepressive_heritable",
  h2_category == "Non-heritable", "nonheritable",
  default = "other"
)]
repr[, intergenic := as.integer(hg38_genes_intergenic == 1)]

for (tis in TISSUES) {
  sub <- repr[tissue == tis]
  cat(sprintf("\n  %s:\n", tis))
  cat(sprintf("    Heritable: %d (repressive: %d, non-repressive: %d)\n",
              sum(sub$h2_category == "Heritable"),
              sum(sub$subgroup == "repressive_heritable"),
              sum(sub$subgroup == "nonrepressive_heritable")))
  cat(sprintf("    Non-heritable: %d\n", sum(sub$h2_category == "Non-heritable")))
  cat(sprintf("    Intergenic heritable: %d (repressive: %d)\n",
              sum(sub$h2_category == "Heritable" & sub$intergenic == 1),
              sum(sub$subgroup == "repressive_heritable" & sub$intergenic == 1)))
  cat(sprintf("    Intergenic non-heritable: %d\n",
              sum(sub$h2_category == "Non-heritable" & sub$intergenic == 1)))
}

## ============================================================
## Step 2: Determine sharing status per VMR
## ============================================================

cat("\n\nDetermining sharing status...\n")

# Load pairwise overlap files (all category) and assign sharing flags
# For each tissue, check if its VMRs appear in the relevant pairwise files.

load_pairwise <- function(tissue1, tissue2) {
  fname <- sprintf("%s_%s_overlap_all.tsv", tissue1, tissue2)
  fpath <- file.path(OVERLAP_DIR, fname)
  dt <- fread(fpath)
  # Harmonize chromosome prefix
  dt[, chromA := paste0("chr", chromA)]
  dt[, chromB := paste0("chr", chromB)]
  dt[, key_a := paste(chromA, startA, endA, sep = ":")]
  dt[, key_b := paste(chromB, startB, endB, sep = ":")]
  dt
}

# Build sharing status for each tissue
# Caudate shares with DLPFC (side A in caudate_dlpfc) and Hippocampus (side A in caudate_hippocampus)
# DLPFC shares with Caudate (side B in caudate_dlpfc) and Hippocampus (side B in hippocampus_dlpfc)
# Hippocampus shares with Caudate (side B in caudate_hippocampus) and DLPFC (side A in hippocampus_dlpfc)

cau_dlp <- load_pairwise("caudate", "dlpfc")
cau_hip <- load_pairwise("caudate", "hippocampus")
hip_dlp <- load_pairwise("hippocampus", "dlpfc")

repr[, shared_with_other1 := FALSE]
repr[, shared_with_other2 := FALSE]

# Caudate
repr[tissue == "Caudate", shared_with_other1 := vmr_key %in% cau_dlp$key_a]
repr[tissue == "Caudate", shared_with_other2 := vmr_key %in% cau_hip$key_a]

# DLPFC
repr[tissue == "DLPFC", shared_with_other1 := vmr_key %in% cau_dlp$key_b]
repr[tissue == "DLPFC", shared_with_other2 := vmr_key %in% hip_dlp$key_b]

# Hippocampus
repr[tissue == "Hippocampus", shared_with_other1 := vmr_key %in% cau_hip$key_b]
repr[tissue == "Hippocampus", shared_with_other2 := vmr_key %in% hip_dlp$key_a]

repr[, n_tissues_shared := 1L + as.integer(shared_with_other1) + as.integer(shared_with_other2)]
repr[, shared_any := as.integer(n_tissues_shared >= 2)]
repr[, shared_all3 := as.integer(n_tissues_shared == 3)]
repr[, tissue_specific := as.integer(n_tissues_shared == 1)]

for (tis in TISSUES) {
  sub <- repr[tissue == tis]
  cat(sprintf("  %s: %d VMRs — tissue-specific: %d (%.0f%%), shared ≥2: %d (%.0f%%), all 3: %d (%.0f%%)\n",
              tis, nrow(sub),
              sum(sub$tissue_specific), 100 * mean(sub$tissue_specific),
              sum(sub$shared_any), 100 * mean(sub$shared_any),
              sum(sub$shared_all3), 100 * mean(sub$shared_all3)))
}

## ============================================================
## Step 3: Compute sharing rates by subgroup
## ============================================================

cat("\nComputing sharing rates by subgroup...\n")

subgroup_levels <- c(
  "repressive_heritable",
  "nonrepressive_heritable",
  "nonheritable",
  "repressive_heritable_intergenic",
  "nonrepressive_heritable_intergenic",
  "nonheritable_intergenic"
)

sharing_rows <- list()

for (tis in TISSUES) {
  sub <- repr[tissue == tis]

  groups <- list(
    repressive_heritable            = sub[subgroup == "repressive_heritable"],
    nonrepressive_heritable         = sub[subgroup == "nonrepressive_heritable"],
    nonheritable                    = sub[h2_category == "Non-heritable"],
    repressive_heritable_intergenic = sub[subgroup == "repressive_heritable" & intergenic == 1],
    nonrepressive_heritable_intergenic = sub[subgroup == "nonrepressive_heritable" & intergenic == 1],
    nonheritable_intergenic         = sub[h2_category == "Non-heritable" & intergenic == 1],
    all_heritable                   = sub[h2_category == "Heritable"],
    all_heritable_intergenic        = sub[h2_category == "Heritable" & intergenic == 1]
  )

  for (grp_name in names(groups)) {
    grp <- groups[[grp_name]]
    n_total <- nrow(grp)
    if (n_total == 0) next
    n_shared <- sum(grp$shared_any)
    n_all3   <- sum(grp$shared_all3)
    ci_shared <- ac_ci(n_shared, n_total)
    ci_all3   <- ac_ci(n_all3, n_total)
    sharing_rows[[length(sharing_rows) + 1]] <- tibble(
      tissue       = tis,
      subgroup     = grp_name,
      n_total      = n_total,
      n_shared     = n_shared,
      sharing_rate = n_shared / n_total,
      ci_lo        = ci_shared$lo,
      ci_hi        = ci_shared$hi,
      n_all3       = n_all3,
      all3_rate    = n_all3 / n_total,
      all3_ci_lo   = ci_all3$lo,
      all3_ci_hi   = ci_all3$hi
    )
  }
}

sharing_df <- bind_rows(sharing_rows)
fwrite(sharing_df, file.path(OUT_DIR, "region_specificity_sharing_rates.tsv"), sep = "\t")
cat("Saved: region_specificity_sharing_rates.tsv\n")

cat("\nSharing rates (shared ≥2 tissues):\n")
sharing_df |>
  dplyr::select(tissue, subgroup, n_total, sharing_rate) |>
  dplyr::mutate(sharing_rate = round(sharing_rate, 3)) |>
  as.data.frame() |>
  print(row.names = FALSE)

## ============================================================
## Steps 4-7: Fisher's exact tests for Q1-Q4
## ============================================================

cat("\nRunning Fisher's exact tests...\n")
fisher_rows <- list()

for (tis in TISSUES) {
  sub <- repr[tissue == tis]

  # Q1: repressive heritable sharing vs all heritable sharing
  rh  <- sub[subgroup == "repressive_heritable"]
  ah  <- sub[h2_category == "Heritable"]
  fisher_rows[[length(fisher_rows) + 1]] <- run_fishers_2x2(
    sum(rh$shared_any), nrow(rh),
    sum(ah$shared_any), nrow(ah),
    tis, "Q1_repressive_vs_all_heritable_sharing",
    "repressive_heritable", "all_heritable"
  )

  # Q2: Caudate-specific vs shared enrichment in repressive marks
  # (run for all tissues for completeness)
  her_specific <- sub[h2_category == "Heritable" & tissue_specific == 1]
  her_shared   <- sub[h2_category == "Heritable" & shared_any == 1]
  if (nrow(her_specific) > 0 & nrow(her_shared) > 0) {
    fisher_rows[[length(fisher_rows) + 1]] <- run_fishers_2x2(
      sum(her_shared$repressive), nrow(her_shared),
      sum(her_specific$repressive), nrow(her_specific),
      tis, "Q2_shared_vs_specific_repressive_enrichment",
      "shared_heritable", "specific_heritable"
    )
    # Also for intergenic
    her_specific_ig <- her_specific[intergenic == 1]
    her_shared_ig   <- her_shared[intergenic == 1]
    if (nrow(her_specific_ig) > 0 & nrow(her_shared_ig) > 0) {
      fisher_rows[[length(fisher_rows) + 1]] <- run_fishers_2x2(
        sum(her_shared_ig$repressive), nrow(her_shared_ig),
        sum(her_specific_ig$repressive), nrow(her_specific_ig),
        tis, "Q2_shared_vs_specific_repressive_intergenic",
        "shared_heritable_intergenic", "specific_heritable_intergenic"
      )
    }
  }

  # Q3: repressive heritable sharing rate vs non-repressive heritable sharing rate
  nrh <- sub[subgroup == "nonrepressive_heritable"]
  fisher_rows[[length(fisher_rows) + 1]] <- run_fishers_2x2(
    sum(rh$shared_any), nrow(rh),
    sum(nrh$shared_any), nrow(nrh),
    tis, "Q3_repressive_vs_nonrepressive_heritable_sharing",
    "repressive_heritable", "nonrepressive_heritable"
  )
  # Q3 intergenic
  rh_ig  <- sub[subgroup == "repressive_heritable" & intergenic == 1]
  nrh_ig <- sub[subgroup == "nonrepressive_heritable" & intergenic == 1]
  if (nrow(rh_ig) > 0 & nrow(nrh_ig) > 0) {
    fisher_rows[[length(fisher_rows) + 1]] <- run_fishers_2x2(
      sum(rh_ig$shared_any), nrow(rh_ig),
      sum(nrh_ig$shared_any), nrow(nrh_ig),
      tis, "Q3_repressive_vs_nonrepressive_heritable_sharing_intergenic",
      "repressive_heritable_ig", "nonrepressive_heritable_ig"
    )
  }

  # Q4: repressive heritable intergenic sharing vs non-heritable intergenic sharing
  nh_ig <- sub[h2_category == "Non-heritable" & intergenic == 1]
  if (nrow(rh_ig) > 0 & nrow(nh_ig) > 0) {
    fisher_rows[[length(fisher_rows) + 1]] <- run_fishers_2x2(
      sum(rh_ig$shared_any), nrow(rh_ig),
      sum(nh_ig$shared_any), nrow(nh_ig),
      tis, "Q4_repressive_heritable_vs_nonheritable_intergenic_sharing",
      "repressive_heritable_ig", "nonheritable_ig"
    )
  }
}

fisher_df <- bind_rows(fisher_rows)
fisher_df$fdr <- p.adjust(fisher_df$p_value, method = "fdr")
fwrite(fisher_df, file.path(OUT_DIR, "region_specificity_fishers.tsv"), sep = "\t")
cat("Saved: region_specificity_fishers.tsv\n")

cat("\nFisher's exact results:\n")
fisher_df |>
  dplyr::mutate(or = round(or, 2), rate_a = round(rate_a, 3), rate_b = round(rate_b, 3),
                fdr = signif(fdr, 2)) |>
  dplyr::select(tissue, question, group_a, rate_a, group_b, rate_b, or, fdr) |>
  as.data.frame() |>
  print(row.names = FALSE)

## ============================================================
## Step 8: Q5 — h² concordance for repressive shared VMRs
## ============================================================

cat("\nComputing h² concordance by repressive status...\n")

# For each tissue pair, load heritable overlap file, merge repressive flags
# from both tissues, and compute Spearman correlation separately for
# repressive-shared and non-repressive-shared pairs.

concordance_rows <- list()

for (pair in TISSUE_PAIRS) {
  t1 <- pair[1]; t2 <- pair[2]
  pair_label <- PAIR_LABELS[paste(t1, t2, sep = "_")]

  # Load heritable pairwise overlaps
  pw <- load_pairwise(t1, t2)

  # Get repressive flags + h2 for tissue 1 and tissue 2
  t1_cap <- unname(TISSUE_LABELS[[t1]])
  t2_cap <- unname(TISSUE_LABELS[[t2]])

  repr_t1 <- repr[tissue == t1_cap, .(vmr_key, h2_unscaled, h2_category,
                                       repressive, intergenic)]
  setnames(repr_t1, c("h2_unscaled", "h2_category", "repressive", "intergenic"),
           c("h2_t1", "h2_cat_t1", "repr_t1", "inter_t1"))

  repr_t2 <- repr[tissue == t2_cap, .(vmr_key, h2_unscaled, h2_category,
                                       repressive, intergenic)]
  setnames(repr_t2, c("h2_unscaled", "h2_category", "repressive", "intergenic"),
           c("h2_t2", "h2_cat_t2", "repr_t2", "inter_t2"))

  # Merge
  pw_merged <- pw |>
    dplyr::mutate(key_a = key_a, key_b = key_b) |>
    dplyr::left_join(as_tibble(repr_t1), by = c("key_a" = "vmr_key")) |>
    dplyr::left_join(as_tibble(repr_t2), by = c("key_b" = "vmr_key")) |>
    dplyr::filter(!is.na(h2_t1), !is.na(h2_t2))

  # Filter to heritable in both tissues
  pw_her <- pw_merged |>
    dplyr::filter(h2_cat_t1 == "Heritable", h2_cat_t2 == "Heritable")

  # Define subsets
  subsets <- list(
    all_heritable       = pw_her,
    repressive_both     = pw_her |> dplyr::filter(repr_t1 == 1, repr_t2 == 1),
    nonrepressive_both  = pw_her |> dplyr::filter(repr_t1 == 0, repr_t2 == 0),
    repressive_either   = pw_her |> dplyr::filter(repr_t1 == 1 | repr_t2 == 1),
    intergenic_both     = pw_her |> dplyr::filter(inter_t1 == 1, inter_t2 == 1),
    repr_intergenic_both = pw_her |> dplyr::filter(repr_t1 == 1, repr_t2 == 1,
                                                    inter_t1 == 1, inter_t2 == 1)
  )

  for (subset_name in names(subsets)) {
    ss <- subsets[[subset_name]]
    n <- nrow(ss)
    if (n < 5) {
      concordance_rows[[length(concordance_rows) + 1]] <- tibble(
        tissue1 = t1_cap, tissue2 = t2_cap, pair_label = pair_label,
        subset = subset_name, n_pairs = n,
        spearman_rho = NA_real_, spearman_p = NA_real_
      )
      next
    }
    sp <- cor.test(ss$h2_t1, ss$h2_t2, method = "spearman", exact = FALSE)
    concordance_rows[[length(concordance_rows) + 1]] <- tibble(
      tissue1 = t1_cap, tissue2 = t2_cap, pair_label = pair_label,
      subset = subset_name, n_pairs = n,
      spearman_rho = sp$estimate, spearman_p = sp$p.value
    )
  }

  # Fisher z-test: repressive_both vs nonrepressive_both
  rho_repr <- subsets$repressive_both |>
    (\(d) if (nrow(d) >= 5) cor(d$h2_t1, d$h2_t2, method = "spearman") else NA_real_)()
  rho_nonr <- subsets$nonrepressive_both |>
    (\(d) if (nrow(d) >= 5) cor(d$h2_t1, d$h2_t2, method = "spearman") else NA_real_)()

  if (!is.na(rho_repr) & !is.na(rho_nonr)) {
    fz <- fisher_z_test(rho_repr, nrow(subsets$repressive_both),
                        rho_nonr, nrow(subsets$nonrepressive_both))
    concordance_rows[[length(concordance_rows) + 1]] <- tibble(
      tissue1 = t1_cap, tissue2 = t2_cap, pair_label = pair_label,
      subset = "z_test_repr_vs_nonrepr", n_pairs = NA_integer_,
      spearman_rho = fz$z_stat, spearman_p = fz$p_value
    )
  }
}

concordance_df <- bind_rows(concordance_rows)
concordance_df$fdr <- p.adjust(concordance_df$spearman_p, method = "fdr")
fwrite(concordance_df, file.path(OUT_DIR, "region_specificity_h2_concordance.tsv"), sep = "\t")
cat("Saved: region_specificity_h2_concordance.tsv\n")

cat("\nh² concordance by repressive status:\n")
concordance_df |>
  dplyr::mutate(spearman_rho = round(spearman_rho, 3), fdr = signif(fdr, 2)) |>
  dplyr::select(pair_label, subset, n_pairs, spearman_rho, fdr) |>
  as.data.frame() |>
  print(row.names = FALSE)

## ============================================================
## Additional: Q2 detailed — caudate-specific repressive enrichment table
## ============================================================

cat("\nQ2 detailed: repressive enrichment by sharing status per tissue\n")

q2_rows <- list()
for (tis in TISSUES) {
  sub <- repr[tissue == tis & h2_category == "Heritable"]
  for (shared_label in c("tissue_specific", "shared_any", "shared_all3")) {
    if (shared_label == "tissue_specific") {
      grp <- sub[tissue_specific == 1]
    } else if (shared_label == "shared_any") {
      grp <- sub[shared_any == 1]
    } else {
      grp <- sub[shared_all3 == 1]
    }
    n <- nrow(grp)
    if (n == 0) next
    n_repr <- sum(grp$repressive)
    ci <- ac_ci(n_repr, n)
    q2_rows[[length(q2_rows) + 1]] <- tibble(
      tissue       = tis,
      sharing      = shared_label,
      n_total      = n,
      n_repressive = n_repr,
      frac         = n_repr / n,
      ci_lo        = ci$lo,
      ci_hi        = ci$hi
    )
  }
}

q2_df <- bind_rows(q2_rows)
fwrite(q2_df, file.path(OUT_DIR, "region_specificity_caudate_enrichment.tsv"), sep = "\t")
cat("Saved: region_specificity_caudate_enrichment.tsv\n")

q2_df |>
  dplyr::mutate(frac = round(frac, 3)) |>
  as.data.frame() |>
  print(row.names = FALSE)

#### Reproducibility ####
cat("\nReproducibility information:\n")
print(Sys.time())
print(proc.time())
options(width = 120)
sessioninfo::session_info()
