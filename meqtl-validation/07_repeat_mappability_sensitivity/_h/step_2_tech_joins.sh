#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --job-name=tech_joins
#SBATCH --output=logs/tech_joins.%j.log

set -euo pipefail
mkdir -p logs
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics
python3 "${ROOT}/meqtl-validation/07_repeat_mappability_sensitivity/_h/04_complete_tech_joins.py" \
  --join-only \
  --min-reciprocal-overlap 0.5
