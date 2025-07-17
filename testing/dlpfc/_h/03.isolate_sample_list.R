# Load required package
library(here)

# Loop through chromosomes 1 to 22
for (chr in 1:22) {
  
  # Construct input and output file paths
  input_file <- here("testing", "dlpfc", "_m", "covs", paste0("chr_", chr), "TOPMed_LIBD.AA.covar")
  output_file <- here("testing", "dlpfc", "_m", "covs", paste0("chr_", chr), "samples.txt")
  
  # Check if input file exists
  if (!file.exists(input_file)) {
    warning("Skipping chromosome ", chr, ": input file not found.")
    next
  }
  
  # Read the .covar file (tab- or space-delimited, with header)
  covar_data <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE)
  
  # Extract the first column (as a vector)
  first_column <- covar_data[[1]]  # Or: covar_data[, 1]
  
  # Write the first column to the output file, one value per line
  writeLines(as.character(first_column), con = output_file)
  
  # Confirmation message
  cat("First column extracted and saved to", output_file, "\n")
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()