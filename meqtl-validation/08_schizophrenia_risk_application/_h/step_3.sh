#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=24G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=scz_cross_prior
#SBATCH --output=logs/scz_cross_prior.%j.log

# Cross-region concordance + ≤5 locus prioritization.
# Submit from meqtl-validation/08_schizophrenia_risk_application/_m/:
#   sbatch ../_h/step_3.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_h"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Cross-region concordance"
python3 "${H}/06_cross_region_concordance.py"

log_message "Locus prioritization (caudate primary)"
python3 "${H}/07_prioritize_loci.py" --primary-region caudate --max-loci 5

log_message "**** Job ends ****"
