#### Identify differentially methylated regions in ###
#### relation to environmental phenotypes ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(data.table)
  library(tidyr)
})

## Function 
get_dmr <- function(pheno_matrix, var) {
  
                                        # Differential DNAm for env phenos
  res <- matrix(NA, nrow=ncol(pheno_matrix), ncol=4)
  colnames(res) <- c("beta","se", "t", "p")
  
  for(i in 1:ncol(pheno_matrix)){
    if(! i %% 100){
      cat(i,"\n")
    }
    d <- as.data.frame(cbind(y=pheno_matrix$meth, pheno_matrix[[var]], 
                             as.factor(pheno_matrix$primarydx), 
                             pheno_matrix$agedeath, 
                             as.factor(pheno_matrix$sex)))
    model = lm(y ~ ., data=d)
    res[i,] = summary(model)$coefficients[2, ]
  }
  res <- as.data.frame(res)
  
                                        # FDR correction
  res$fdr <- p.adjust(res$p, method="fdr")
  
  res <- as.data.frame(cbind(pheno_matrix, res))
  colnames(res)[1:3] <- c("chr", "start", "end")
  
                                        # Output significant dmrs
  out_dmr <- file.path(out_path, paste0(var, "_dmr.csv"))
  write.csv(res, out_dmr)

  return(res)
}

## Main
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

#meth_matrix <- pheno_matrix %>%
#  select("brnum", "feature_id", "meth") %>%
#  pivot_wider(names_from = "feature_id", values_from = "meth")

for (env in vars_to_include) {
  res <- get_dmr(pheno_matrix, env)
  
                                        # Count significant dmrs
  print(paste("At FDR < 0.05 there are", sum(res$fdr < 0.05), "signficant DMRs"))
  print(paste("At FDR < 0.1 there are", sum(res$fdr < 0.1), "signficant DMRs"))
  print(paste("At p < 0.05 there are", sum(res$p < 0.05), "signficant DMRs"))
  
}

#### Reproducibility ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()