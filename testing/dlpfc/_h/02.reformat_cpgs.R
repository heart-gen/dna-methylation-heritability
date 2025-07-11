# Load the 'here' package
library(here)

# Input and output file paths using 'here'
input_file <- here("testing/dlpfc/_m/cpg/chr_1/cpg_pos.txt")
output_file <- here("testing/dlpfc/_m/cpg/chr_1/reformatted_cpg_pos.txt")

# Read input numbers
numbers <- as.numeric(readLines(input_file))

# Modify each line
modified_lines <- mapply(function(x, i) {
  paste0("22\t", x - 500000, "\t", x + 500000, "\tR", i)
}, numbers, seq_along(numbers))

# Write to output file
writeLines(modified_lines, output_file)

cat("Modified file written to", output_file, "\n")