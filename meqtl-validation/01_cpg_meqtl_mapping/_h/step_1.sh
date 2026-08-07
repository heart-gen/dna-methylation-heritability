#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_cpg_preflight
#SBATCH --output=logs/cpg_preflight.%A_%a.log

# CPU: sample inclusion + covariates
# Submit from meqtl-validation/01_cpg_meqtl_mapping/_m/:
#   mkdir -p logs && sbatch ../_h/step_1.sh
# Single region: sbatch --export=ALL,REGION=caudate --array=0 ../_h/step_1.sh

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
echo "**** QUEST info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Node name: ${SLURM_NODENAME:-N/A}"
echo "Array task: ${SLURM_ARRAY_TASK_ID:-N/A}"

module purge
module list
mkdir -p logs

REGIONS=(caudate dlpfc hippocampus)
if [[ -z "${REGION:-}" ]]; then
  REGION="${REGIONS[${SLURM_ARRAY_TASK_ID}]}"
fi
POPULATION="${POPULATION:-AA}"

log_message "Preflight + covariates for ${REGION} (${POPULATION})"
python3 "../_h/00_preflight.py" --region "${REGION}" --population "${POPULATION}"
python3 "../_h/01_prepare_covariates.py" --region "${REGION}" --population "${POPULATION}"

log_message "**** Job ends ****"
