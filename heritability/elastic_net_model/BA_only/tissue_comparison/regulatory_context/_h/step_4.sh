#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=regctx_prox
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/regctx_prox.%A_%a.log
#SBATCH --array=0-2

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

log_message "Calculating feature proximity"

RUNS=("caudate AA" "dlpfc AA" "hippocampus AA")
read -r TISSUE POP <<< "${RUNS[$SLURM_ARRAY_TASK_ID]}"
Rscript "../_h/04.feature_proximity.R" \
  BA_only "$TISSUE" "$POP"

conda deactivate
log_message "**** Job ends ****"
