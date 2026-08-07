#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --array=0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_vmr_models
#SBATCH --output=logs/vmr_models.%A_%a.log

# Submit from meqtl-validation/02_vmr_meqtl_burden/_m/ after step_1:
#   sbatch ../_h/step_2.sh

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
mkdir -p logs

REGIONS=(caudate dlpfc hippocampus)
if [[ -z "${REGION:-}" ]]; then
  REGION="${REGIONS[${SLURM_ARRAY_TASK_ID}]}"
fi
POPULATION="${POPULATION:-AA}"

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h"
if [[ "${POPULATION}" == "AA" ]]; then
  BURDEN="${ROOT}/meqtl-validation/02_vmr_meqtl_burden/_m/${REGION}/vmr_meqtl_burden.tsv.gz"
else
  BURDEN="${ROOT}/meqtl-validation/02_vmr_meqtl_burden/_m/${POPULATION}/${REGION}/vmr_meqtl_burden.tsv.gz"
fi

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Fitting VMR burden models for ${REGION} (${POPULATION})"
python3 "${H}/02_fit_burden_models.py" --region "${REGION}" --burden-tsv "${BURDEN}"

log_message "**** Job ends ****"
