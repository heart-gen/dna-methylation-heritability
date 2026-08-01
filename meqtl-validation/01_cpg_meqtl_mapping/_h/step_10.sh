#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_covsens_cmp
#SBATCH --output=logs/covsens_compare.%j.log

# Compare finished covariate models and write lock decision.
# Submit from _m/ after step_9 jobs complete:
#   sbatch --export=ALL,REGION=caudate ../_h/step_10_compare_models.sh
# Optional: build M4 after choosing k, then re-run step_9 for M4 only.

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

REGION="${REGION:-caudate}"
PEER_K_FOR_M4="${PEER_K_FOR_M4:-0}"
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

if [[ "${PEER_K_FOR_M4}" != "0" ]]; then
  log_message "Building M4 with k=${PEER_K_FOR_M4} for ${REGION}"
  python3 "${H}/10_prepare_covariate_models.py" \
    --region "${REGION}" \
    --peer-k-for-m4 "${PEER_K_FOR_M4}"
fi

log_message "Comparing covariate models for ${REGION}"
python3 "${H}/11_compare_covariate_models.py" --region "${REGION}"

log_message "Writing lock decision (pilot-driven)"
python3 "${H}/12_lock_decision.py" --pilot-region "${REGION}"

log_message "**** Job ends ****"
