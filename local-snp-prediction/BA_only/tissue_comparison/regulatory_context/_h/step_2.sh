#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=regctx_arch
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/regctx_arch.%A_%a.log
#SBATCH --array=0-5

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

echo "**** QUEST info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

module purge
module list

mkdir -p ../_m/logs
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

RUNS=(
  "caudate expression AA 250000"
  "caudate psi AA 250000"
  "dlpfc expression AA 250000"
  "dlpfc psi AA 250000"
  "hippocampus expression AA 250000"
  "hippocampus psi AA 250000"
)

log_message "Comparing architecture"

read -r TISSUE MODALITY POP WINDOW <<< "${RUNS[$SLURM_ARRAY_TASK_ID]}"
Rscript "../_h/02.compare_architecture.R" \
  BA_only "$TISSUE" "$MODALITY" "$POP" "$WINDOW"

conda deactivate
log_message "**** Job ends ****"
