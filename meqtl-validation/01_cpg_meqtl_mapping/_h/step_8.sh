#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_latent
#SBATCH --output=logs/covsens_latent.%A_%a.log

# Estimate methylation latent factors + build covariate matrices.
# Submit from meqtl-validation/01_cpg_meqtl_mapping/_m/:
#   mkdir -p logs && sbatch ../_h/step_8_latent.sh
# Single region:
#   sbatch --export=ALL,REGION=caudate --array=0 ../_h/step_8_latent.sh

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

REGIONS=(caudate dlpfc hippocampus)
if [[ -z "${REGION:-}" ]]; then
  REGION="${REGIONS[${SLURM_ARRAY_TASK_ID}]}"
fi

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Estimating latent factors for ${REGION}"
python3 "${H}/09_estimate_latent_factors.py" --region "${REGION}"

log_message "Building covariate model matrices for ${REGION}"
python3 "${H}/10_prepare_covariate_models.py" --region "${REGION}"

log_message "**** Job ends ****"
