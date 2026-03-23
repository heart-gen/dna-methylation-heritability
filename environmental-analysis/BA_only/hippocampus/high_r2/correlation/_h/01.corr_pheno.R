#### Associate environmental phenotypes with low-heritability VMRs ####

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(tidyr)
})

## Function           
filter_sites <- function(enet) {
  vmr <- na.omit(enet)
  vmr <- vmr %>% 
    dplyr::select(chrom, start, end, h2_unscaled, r_squared_cv) %>% 
    mutate(h2_category = case_when(
      r_squared_cv <= 0.75 ~ "Low prediction",
      h2_unscaled < 0.1 & r_squared_cv > 0.75 ~ "Non-heritable",
      h2_unscaled >= 0.1 & r_squared_cv > 0.75 ~ "Heritable"
    ),
    h2_category = factor(h2_category,
                         levels = c("Heritable", "Non-heritable", "Low prediction"))
    )
  
  return(vmr)
}

clean_pheno <- function(pheno_file_path, tissue, vars){
  pheno_df <- fread(pheno_file_path, header = TRUE) |>
    dplyr::select(all_of(vars_to_include)) |>
    filter(agedeath > 17, region == tissue) |>
    mutate(education = case_when(
      education %in% c("7th", "8th", "Less than 7th",
                       "9th", "10th", "11th", "12th") ~ "less_than_hs",
      education %in% c("H.S. diploma","GED") ~ "hs",
      education %in% c("1 yr college", "3 yrs college", "Associate's or 2 yrs college",
                       "Bachelor's", "Master's", "JD", "PhD") ~ "more_than_hs"
    )
    ) |>
    mutate(marital_status = case_when(
      marital_status %in% c("Single") ~ "single",
      marital_status %in% c("Married") ~ "married",
      marital_status %in% c("Divorced", "Separated", "Widowed") ~ "previously_married"
    )
    ) |>
    rename(age = agedeath, dx = primarydx) |>
    mutate_if(is.character, as.factor)
  
  return(pheno_df)
}

merge_meth <- function(meth_files){
  meth_list <- vector("list", length(meth_files))
  for (i in seq_along(meth_files)) {
    file  <- meth_files[i]
                                        # Get pos from filename
    chr   <- basename(dirname(file))
    pos   <- strsplit(sub("_meth\\.phen$", "", basename(file)), "_")[[1]]
    df    <- fread(meth_files[i], select = c("V1", "V3"))
    colnames(df) <- c("brnum", "meth")
    
                                        # Add in vmr pos
    df <- df %>%
      mutate(chr   = as.integer(sub("chr_","", chr)),
             start = as.integer(pos[1]),
             end   = as.integer(pos[2]),
             feature_id = paste0("VMR", i))
    meth_list[[i]] <- df
  }
                                        # Bind meth matrix
  meth_df <- rbindlist(meth_list, use.names = TRUE, fill = TRUE)
  
  return(meth_df)
}

save_plot <- function(p, fn, w=8, h=6) {
  for (ext in c(".png", ".pdf")) {
    ggsave(filename = paste0(fn, ext), plot = p, width = w, height = h)
  }
}

## Main
tissue <- c("hippocampus")

out_path <- here("environmental-analysis", "BA_only", 
                 "hippocampus", "high_r2", "correlation", "_m")
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

                                        # Read in summary table
enet_file <- here("heritability/elastic_net_model/BA_only/", 
                  paste0(tissue, "/_m/", tissue, "_summary_elastic-net.tsv"))
enet <- read.table(enet_file, sep = "\t", header = TRUE)

vmr <- filter_sites(enet)

                                        # Define variables of interest
vars_to_include <- c(
  "brnum", "agedeath", "sex","primarydx","region","smoking","codeine","morphine",
  "cocaine","ethanol","antipsychotics","nicotine","amphetamines",
  "education","marital_status","hx_sexual_abuse","hx_physical_abuse",
  "hx_other_trauma","hx_military_service","fsiq"
)

                                        # Get phenotype data
pheno_file_path <- here("inputs/phenotypes/_m/phenotypes-AA.tsv")
pheno <- clean_pheno(pheno_file_path, tissue, vars_to_include)
                    
                                        # Add global ances
f_ances   <- here("inputs", "genetic-ancestry",
                  "structure.out_ancestry_proportion_raceDemo_compare")
ances     <- fread(f_ances) %>% filter(group == "AA")
pheno <- left_join(pheno, ances, by = c("brnum" = "id")) %>%
  rename(afr_ances = Afr)

                                        # Get meth matrix
meth_file_path <- here("vmr-analysis/hippocampus/_m/vmr")
meth_files     <- list.files(path = meth_file_path, pattern = "_meth\\.phen$", 
                         recursive = TRUE, full.names = TRUE)
meth_df        <- merge_meth(meth_files)
meth_df <- meth_df %>% mutate(chr = as.character(chr))

                                        # Merge with h2 groups and 
                                        # environmental variables
groups <- meth_df |>
  inner_join(vmr, by = c("chr" = "chrom", "start", "end"))
merged <- groups |>
  inner_join(pheno, "brnum") |>
  arrange(feature_id)

                                        # Write df to file
out_table <- file.path(out_path, paste0("vmr_env_assoc-AA.tsv.gz"))
fwrite(merged, out_table, sep = "\t", quote = FALSE)


####### Covariate testing #########

                                        # Define covariates
covars <- "age + sex + dx + afr_ances"

                                        # Initialize results matrix
cov_results <- tibble() 
feature_ids <- unique(merged$feature_id)

for (vmr in feature_ids) {
  print(paste("Running null logistic regression on", vmr))
  
  vmr_merged <- merged %>% filter(feature_id == vmr)
  
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
pos <- merged %>%
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
testing_envs <- c(
  "smoking", "codeine", "morphine", "cocaine", "ethanol", "antipsychotics", 
  "nicotine", "amphetamines", "education", "marital_status", "hx_sexual_abuse",
  "hx_physical_abuse", "hx_other_trauma", "hx_military_service", "fsiq"
)
  
for (env in testing_envs) {
  merged %>% 
    group_by(across(all_of(env)), h2_category) %>%
    summarize(
      count = n(),
      mean = mean(meth),
      sd = sd(meth)
    )
  
  merged_env <- merged %>% drop_na(env)
  
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
      mutate(beta = Estimate, se = `Std. Error`,
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