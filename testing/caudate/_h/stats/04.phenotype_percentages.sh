#!/bin/bash
#SBATCH -A p32505
#SBATCH -p short
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --job-name=phenotype_percentages
#SBATCH --output=logs/phenotype_percentages.out
#SBATCH --error=logs/phenotype_percentages.err

python /projects/p32505/users/elisa/testing/caudate/_h/phenotype_percentages.py