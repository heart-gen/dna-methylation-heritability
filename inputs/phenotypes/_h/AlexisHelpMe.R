##############

install.packages("dyplr")
install.packages("tidyverse")
install.packages("here")

library(dplyr)
library(tidyverse)
library(here)

phenotype_data <- read.table(here("inputs/phenotypes/_m/phenotypes-AA.tsv"), 
                             sep = "\t", header = TRUE)

phenotype_data <- phenotype_data %>%
  filter(agedeath >= 17)

phenotypes_summary <- phenotype_data %>%
  group_by(region) %>%
  summarize(
    patients = n(), 
    male = sum(sex == "M"), 
    female = sum(sex == "F"), 
    
    SCZ = sum(primarydx == "Schizo"),
    NT = sum(primarydx == "Control"),
    average_age = mean("agedeath"),
    
  )
print(phenotypes_summary)

output_path <- here("inputs", "phenotypes", "_m")

if(!dir.exist(output_path)) {
  dir.create(output_path, recursive = TRUE) }


write.csv(phenotypes_summary, file = file.path(output_path, "demographics_summary.csv"))


################### REPRO

print("Reproducibility information")
sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

