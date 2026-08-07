#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_cpg_genotypes
#SBATCH --output=logs/cpg_genotypes.%A_%a.log

# CPU: filter genotypes to analysis samples (plink2)
# Requires step_1 inclusion list. Submit from _m/:
#   sbatch ../_h/step_3.sh

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

log_message "**** Loading conda environment ****"
# shellcheck disable=SC1091
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics
command -v plink2 >/dev/null || { log_message "ERROR: plink2 missing after activating genomics env"; exit 1; }

REGIONS=(caudate dlpfc hippocampus)
if [[ -z "${REGION:-}" ]]; then
  REGION="${REGIONS[${SLURM_ARRAY_TASK_ID}]}"
fi
POPULATION="${POPULATION:-AA}"

log_message "Filtering genotypes for ${REGION} (${POPULATION})"
bash "../_h/03_prepare_genotypes.sh" "${REGION}" "${POPULATION}"

conda deactivate
log_message "**** Job ends ****"
