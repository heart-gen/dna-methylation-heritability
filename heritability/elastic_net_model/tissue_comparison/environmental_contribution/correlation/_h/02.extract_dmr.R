#### Identify differentially methylated regions in ###
#### relation to environmental phenotypes ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(limma)
})

## Function 
get_dmr <- function(pheno_matrix, var) {
  
                                        # Format matrices
  meth_levels <- pheno_matrix %>% drop_na(var) %>%
    select("brnum", "feature_id", "meth") %>%
    pivot_wider(names_from = "feature_id", values_from = "meth")
  
  meth_levels <- meth_levels[ , -1]
  
  pheno_matrix <- pheno_matrix %>% drop_na(var) %>%
    distinct(brnum, .keep_all = TRUE)
                                        # Differential DNAm for env phenos
  meth_t <- t(meth_levels)
  
  design <- cbind(1, pheno_matrix[[var]], 
                  as.factor(pheno_matrix$primarydx), 
                  pheno_matrix$agedeath, 
                  as.factor(pheno_matrix$sex))
  fit    <- limma::lmFit(meth_t, design)
  fit    <- eBayes(fit)
  
  
  res <- as.data.frame(fit)
  
                                        # FDR correction
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
pheno_matrix_fn <- here("heritability", "elastic_net_model", "tissue_comparison",
                        "environmental_contribution", "correlation", "_m",
                        "vmr_env_assoc-AA.tsv.gz")
pheno_matrix <- fread(pheno_matrix_fn)
vmr_ids <- pheno_matrix %>% distinct(feature_id)

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
  
  out_dmr <- file.path(out_path, paste0(env, "_dmr.csv"))
  write.csv(merged, out_dmr)
}

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()