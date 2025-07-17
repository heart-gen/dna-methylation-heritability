# Load required package
library(here)

# Loop through chromosomes 1 to 22
for (chr in 1:22) {
  
  # Construct input and output file paths
  input_file <- here("testing", "dlpfc", "_m", "cpg", paste0("chr_", chr), "cpg_pos.txt")
  output_file <- here("testing", "dlpfc", "_m", "cpg", paste0("chr_", chr), "cpg_pos_reformatted.txt")
  
  # Check if input file exists
  if (!file.exists(input_file)) {
    warning("Skipping chromosome ", chr, ": input file not found.")
    next
  }
  
  # Read all lines and remove the first 2
  all_lines <- readLines(input_file)
  if (length(all_lines) <= 2) {
    warning("Skipping chromosome ", chr, ": file has less than or equal to 2 lines.")
    next
  }
  filtered_lines <- all_lines[-c(1,2)]
  
  # Convert to numeric
  numbers <- as.numeric(filtered_lines)
  
  # Modify each line
  modified_lines <- mapply(function(x, i) {
    paste0(chr, "\t", x - 500000, "\t", x + 500000, "\tR", i)
  }, numbers, seq_along(numbers))
  
  # Write output
  writeLines(modified_lines, output_file)
  cat("Wrote chromosome", chr, "to", output_file, "\n")
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()