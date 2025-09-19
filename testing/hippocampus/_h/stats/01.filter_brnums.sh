#!/bin/bash
#SBATCH -A p32505
#SBATCH -p short
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --job-name=filter_brnums
#SBATCH --output=logs/filter_brnums.out
#SBATCH --error=logs/filter_brnums.err

python /projects/p32505/users/elisa/testing/hippocampus/_h/01.filter_brnums.py