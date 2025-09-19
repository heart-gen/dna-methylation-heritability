#!/bin/bash
#SBATCH -A p32505
#SBATCH -p short
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --job-name=count_sexes
#SBATCH --output=logs/count_sexes.out
#SBATCH --error=logs/count_sexes.err

python /projects/p32505/users/elisa/testing/dlpfc/_h/count_sexes.py