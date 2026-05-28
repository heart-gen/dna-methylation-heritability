#### Associate environmental phenotypes with low-heritability VMRs ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

source(here("environmental-analysis", "BA_only", "tissue_compare", "_h", "discovery_exposures.R"))
N_QUINTILES <- 5

## Function
create_quintiles <- function(df) {
  # Remove low prediction VMRs
  df <- df %>%
    filter(h2_category %in% c("Heritable", "Non-heritable"))
  
  breaks <- quantile(df$h2_unscaled,
                     probs = seq(0, 1, 1 / N_QUINTILES),
                     na.rm = TRUE) |>
    unique()
  n_bins <- length(breaks) - 1

  if (n_bins < 1) {
    warning(sprintf(
      "Insufficient unique h2 values to compute quantile bins"
    ))
    return(tibble())
  }

  if (n_bins < N_QUINTILES) {
    warning(sprintf(
      "Reduced h2 bins from %d to %d because quantile breaks were not unique",
      N_QUINTILES, n_bins
    ))
  }

  df |>
    mutate(h2_quintile = cut(h2_unscaled,
                             breaks = breaks,
                             labels = paste0("Q", seq_len(n_bins)),
                             include.lowest = TRUE)) |>
    filter(!is.na(h2_quintile))
}

save_plot <- function(p, fn, w=8, h=6) {
  for (ext in c(".png", ".pdf")) {
    ggsave(filename = paste0(fn, ext), plot = p, width = w, height = h)
  }
}

## Main
tissue <- c("caudate")

out_path <- here("environmental-analysis", "BA_only", "quintile_analysis"
                 paste0(tissue, "/correlation/_m"))
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

                                        # Get phenotype matrix
pheno_matrix_fn <- here("environmental-analysis", "BA_only", paste0(tissue, 
                        "/correlation/_m/vmr_env_assoc-AA.tsv.gz"))
pheno_matrix <- fread(pheno_matrix_fn, na.strings = c(NA, ""))
report_env_var_coverage(pheno_matrix, get_discovery_env_vars())
pheno_matrix_split <- create_quintiles(pheno_matrix)

####### Covariate testing #########

                                        # Define covariates
covars <- "age + sex + dx + afr_ances"

                                        # Get unique quintiles
valid_quintiles <- sort(unique(pheno_matrix_split$h2_quintile))

for (quintile in valid_quintiles) {
  print(paste("Running regression models in", quintile))
  
  subdir <- paste(quintile)
  
  # create output directories if they  
  # don't exist
  subdir_path <- file.path(out_path, subdir)
  if (!dir.exists(subdir_path)) {
    dir.create(subdir_path, recursive = TRUE)
  }

                                          # Initialize results matrix
  cov_results <- tibble()
  pheno_quintile <- pheno_matrix_split %>% filter(h2_quintile == quintile)
  feature_ids <- unique(pheno_quintile$feature_id)

  for (vmr in feature_ids) {
    print(paste("Running null linear regression on", vmr))

    vmr_merged <- pheno_quintile %>% filter(feature_id == vmr)

    # Test control variables
    cov_model <- as.formula(paste("meth ~", covars))
    cov_lm <- lm(cov_model, data = vmr_merged)

    cov_summary <- as.data.frame(summary(cov_lm)$coefficients)
    cov_summary$var <- rownames(cov_summary)
    cov_summary <- cov_summary %>%
      filter(var != "(Intercept)") %>%
      mutate(feature_id = vmr) %>%
      rename(beta = Estimate, se = `Std. Error`,
            t = `t value`, p = `Pr(>|t|)`)

    cov_results <- bind_rows(cov_results, cov_summary)
  }

                                          # FDR correction
  cov_results <- cov_results %>% as.data.frame() %>%
    mutate(fdr = p.adjust(p, method = "fdr"))

                                          # Add positions back in
  pos <- pheno_quintile %>%
    dplyr::select(feature_id, chr, start, end) %>%
    distinct()
  cov_results <- pos %>%
    inner_join(cov_results, "feature_id") %>%
    mutate(var = recode(var, "sexM" = "sex"),
          var = recode(var, "dxSchizo" = "dx"))

  cov_names <- unique(cov_results$var)

  for (cov in cov_names){
    cov_filtered <- cov_results %>% filter(startsWith(var, cov))

    # Write results
    out_cov_lm <- file.path(subdir_path, paste0(cov, "_linear.csv.gz"))
    fwrite(cov_filtered, out_cov_lm)
  }

  ####### Environmental testing #########

                                          # Define environmental vars
  testing_envs <- get_discovery_env_vars()

  for (env in testing_envs) {
    pheno_quintile %>%
      group_by(across(all_of(env)), h2_category) %>%
      summarize(
        count = n(),
        mean = mean(meth),
        sd = sd(meth)
      )

    merged_env <- pheno_quintile %>% drop_na(all_of(env))

    # Initialize results matrix
    env_results <- tibble()

    # Grouped logistic regression
    for (vmr in feature_ids) {
      print(paste("Running linear regression on", vmr, "as a function of", env))

      vmr_merged <- merged_env %>% filter(feature_id == vmr)

      # Test environmental variables
      env_model <- as.formula(paste("meth ~", env, "+", covars))
      env_lm <- lm(env_model, data = vmr_merged)

      env_summary <- as.data.frame(summary(env_lm)$coefficients)
      env_summary$var <- rownames(env_summary)
      env_summary <- env_summary %>%
        filter(grepl(paste0("^", env), var)) %>%
        mutate(feature_id = vmr, env = env) %>%
        rename(beta = Estimate, se = `Std. Error`,
              t = `t value`, p = `Pr(>|t|)`)

      env_results <- bind_rows(env_results, env_summary)
    }
    # FDR correction
    env_results <- env_results %>% as.data.frame() %>%
      mutate(fdr = p.adjust(p, method = "fdr"))

    # Count significant VMRs
    print(paste("At FDR < 0.05 there are", sum(env_results$fdr < 0.05), "signficant VMRs"))
    print(paste("At FDR < 0.1 there are", sum(env_results$fdr < 0.1), "signficant VMRs"))
    print(paste("At p < 0.05 there are", sum(env_results$p < 0.05), "signficant VMRs"))

    # Add positions back in
    pos <- merged_env %>%
      dplyr::select(feature_id, chr, start, end) %>%
      distinct()
    env_results <- pos %>%
      inner_join(env_results, "feature_id")

    # Write results
    out_lm <- file.path(subdir_path, paste0(env, "_linear.csv.gz"))
    fwrite(env_results, out_lm)
  }
}

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
