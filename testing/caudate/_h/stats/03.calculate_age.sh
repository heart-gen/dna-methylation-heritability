#!/bin/bash
#SBATCH -A p32505
#SBATCH -p short
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --job-name=calculate_age
#SBATCH --output=logs/calculate_age.out
#SBATCH --error=logs/calculate_age.err

python /projects/p32505/users/elisa/testing/caudate/_h/calculate_age.py