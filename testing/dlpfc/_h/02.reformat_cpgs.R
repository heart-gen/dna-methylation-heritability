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
  
  # Read input numbers
  numbers <- as.numeric(readLines(input_file))
  
  # Modify each line
  modified_lines <- mapply(function(x, i) {
    paste0(chr, "\t", x - 500000, "\t", x + 500000, "\tR", i)
  }, numbers, seq_along(numbers))
  
  # Write output
  writeLines(modified_lines, output_file)
  cat("Wrote chromosome", chr, "to", output_file, "\n")
}
