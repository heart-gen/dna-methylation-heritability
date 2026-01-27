#### Identify differentially methylated regions in ###
#### relation to environmental phenotypes ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(limma)
  library(tibble)
})

## Function 
get_null_dmr <- function(pheno_matrix) {
  
  # Format matrices
  meth_levels <- pheno_matrix %>%
    select("brnum", "feature_id", "meth") %>%
    pivot_wider(names_from = "feature_id", values_from = "meth")
  
  meth_levels <- meth_levels[ , -1]
  
  pheno_matrix <- pheno_matrix %>%
    distinct(brnum, .keep_all = TRUE) %>%
    mutate(education = as.numeric(as.factor(education))) %>%
    mutate(marital_status = as.numeric(as.factor(marital_status)))
  
  # Differential DNAm for env phenos
  meth_t <- t(meth_levels)
  
  design <- cbind(1, as.factor(pheno_matrix$primarydx), 
                  pheno_matrix$agedeath, 
                  as.factor(pheno_matrix$sex))
  fit    <- limma::lmFit(meth_t, design)
  fit    <- eBayes(fit)
  
  coef_beta <- as.data.frame(fit$coefficients) %>% 
    rownames_to_column("feature_id") %>%
    pivot_longer(cols = -feature_id, names_to = "var", values_to = "coefficients")
  coef_p <- as.data.frame(fit$p.value) %>% 
    rownames_to_column("feature_id") %>%
    pivot_longer(cols = -feature_id, names_to = "var", values_to = "p.value")
  
  coef_df <- left_join(coef_beta, coef_p, by = c("feature_id", "var"))
  
  # FDR correction
  coef_res <- coef_df %>%
    group_by(var) %>%
    mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
    ungroup()
  
  return(coef_res)
}

get_dmr <- function(pheno_matrix, var) {
  
                                        # Format matrices
  meth_levels <- pheno_matrix %>% drop_na(var) %>%
    select("brnum", "feature_id", "meth") %>%
    pivot_wider(names_from = "feature_id", values_from = "meth")
  
  meth_levels <- meth_levels[ , -1]
  
  pheno_matrix <- pheno_matrix %>% drop_na(var) %>%
    distinct(brnum, .keep_all = TRUE) %>%
    mutate(education = as.numeric(as.factor(education))) %>%
    mutate(marital_status = as.numeric(as.factor(marital_status)))
  
                                        # Differential DNAm for env phenos
  meth_t <- t(meth_levels)
  
  design <- cbind(1, pheno_matrix[[var]], 
                  as.factor(pheno_matrix$primarydx), 
                  pheno_matrix$agedeath, 
                  as.factor(pheno_matrix$sex))
  fit    <- limma::lmFit(meth_t, design)
  fit    <- eBayes(fit)
  
                                        # FDR correction
  res <- as.data.frame(fit)
  res$fdr <- p.adjust(res$p.value.x2, method="fdr")

  return(res)
}

## Main
                                        # Create output path
out_path <- here("heritability", "elastic_net_model", "tissue_comparison",
                 "environmental_contribution", "correlation", "_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

                                        # Define variables of interest
vars_to_include <- c(
  "smoking","codeine","morphine", "cocaine","ethanol","antipsychotics",
  "nicotine","amphetamines", "education","marital_status","hx_sexual_abuse",
  "hx_physical_abuse", "hx_other_trauma","hx_military_service","fsiq"
)

                                        # Get phenotype matrix
pheno_matrix_fn <- here("environmental-analysis", "correlation", "_m",
                        "vmr_env_assoc-AA.tsv.gz")
pheno_matrix <- fread(pheno_matrix_fn, na.strings = c(NA, ""))

                                        # Add global ances
f_ances   <- here("inputs", "genetic-ancestry",
                  "structure.out_ancestry_proportion_raceDemo_compare")
ances     <- fread(f_ances) %>% filter(group == "AA")
pheno_matrix <- left_join(pheno_matrix, ances, by = c("brnum" = "id"))

                                        # Get VMR IDs
vmr_ids <- pheno_matrix %>% distinct(feature_id)

print(paste("Extracting DMRs related to covariates"))

coef_res <- get_null_dmr(pheno_matrix)

# Count significant dmrs
summary <- coef_res %>%
  group_by(var) %>%
  summarise(sig_fdr = sum(fdr < 0.05),
            sig_p = sum(p.value < 0.05))
print(summary)

# Output significant dmrs
cov_merged <- as.data.frame(cbind(vmr_ids, coef_res)) %>%
  left_join(select(pheno_matrix, chr, start, end, feature_id), multiple = "first") %>%
  select("feature_id", "chr", "start", "end", everything())

out_cov_dmr <- file.path(out_path, "cov_dmr.csv.gz")
fwrite(cov_merged, out_cov_dmr)

# Run DMR analysis for environmental vars
for (env in vars_to_include) {
  
  print(paste("Extracting DMRs related to:", env))
  
  res <- get_dmr(pheno_matrix, env)
  
                                        # Count significant dmrs
  print(paste("At FDR < 0.05 there are", sum(res$fdr < 0.05), "signficant DMRs"))
  print(paste("At FDR < 0.1 there are", sum(res$fdr < 0.1), "signficant DMRs"))
  print(paste("At p < 0.05 there are", sum(res$p.value.x2 < 0.05), "signficant DMRs"))
  
                                        # Output significant dmrs
  merged <- as.data.frame(cbind(vmr_ids, res)) %>%
    left_join(select(pheno_matrix, chr, start, end, feature_id), multiple = "first") %>%
    select("feature_id", "chr", "start", "end", everything())
  
  out_dmr <- file.path(out_path, paste0(env, "_dmr.csv.gz"))
  fwrite(merged, out_dmr)
}

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()