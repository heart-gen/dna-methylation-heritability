#### Associate environmental phenotypes with low-heritability VMRs ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(data.table)
})

source(here("environmental-analysis", "BA_only", "tissue_compare", "_h", "discovery_exposures.R"))

## Function
save_plot <- function(p, fn, w=8, h=6) {
  for (ext in c(".png", ".pdf")) {
    ggsave(filename = paste0(fn, ext), plot = p, width = w, height = h)
  }
}

## Main
tissue <- c("dlpfc")

out_path <- here("environmental-analysis", "BA_only", 
                 paste0(tissue, "/correlation/_m"))
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

                                        # Get phenotype matrix
pheno_matrix_fn <- here("environmental-analysis", "BA_only", paste0(tissue, 
                        "/correlation/_m/vmr_env_assoc-AA.tsv.gz"))
pheno_matrix <- fread(pheno_matrix_fn, na.strings = c(NA, ""))
report_env_var_coverage(pheno_matrix, get_discovery_env_vars())

####### Covariate testing #########

                                        # Define covariates
covars <- "age + sex + dx + afr_ances"

                                        # Initialize results matrix
cov_results <- tibble()
feature_ids <- unique(pheno_matrix$feature_id)

for (vmr in feature_ids) {
  print(paste("Running null logistic regression on", vmr))

  vmr_merged <- pheno_matrix %>% filter(feature_id == vmr)

  # Test control variables
  cov_model <- as.formula(paste("meth ~", covars))
  cov_logit <- glm(cov_model, data = vmr_merged)

  cov_summary <- as.data.frame(summary(cov_logit)$coefficients)
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
pos <- pheno_matrix %>%
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
  out_cov_logit <- file.path(out_path, paste0(cov, "_logit.csv.gz"))
  fwrite(cov_filtered, out_cov_logit)
}

####### Environmental testing #########

                                        # Define environmental vars
testing_envs <- get_discovery_env_vars()

for (env in testing_envs) {
  pheno_matrix %>%
    group_by(across(all_of(env)), h2_category) %>%
    summarize(
      count = n(),
      mean = mean(meth),
      sd = sd(meth)
    )

  merged_env <- pheno_matrix %>% drop_na(env)

  # Initialize results matrix
  env_results <- tibble()

  # Grouped logistic regression
  for (vmr in feature_ids) {
    print(paste("Running logistic regression on", vmr, "as a function of", env))

    vmr_merged <- merged_env %>% filter(feature_id == vmr)

    # Test environmental variables
    env_model <- as.formula(paste("meth ~", env, "+", covars))
    env_logit <- glm(env_model, data = vmr_merged)

    env_summary <- as.data.frame(summary(env_logit)$coefficients)
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
  out_logit <- file.path(out_path, paste0(env, "_logit.csv.gz"))
  fwrite(env_results, out_logit)
}

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()