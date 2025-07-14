#!/bin/bash
#SBATCH --job-name=modify_lines
#SBATCH --output=modify_lines.out
#SBATCH --error=modify_lines.err
#SBATCH --time=01:00:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task=1

# Load R module (adjust as needed)
module load R

# Run the R script
Rscript 02.reformat_cpgs.R