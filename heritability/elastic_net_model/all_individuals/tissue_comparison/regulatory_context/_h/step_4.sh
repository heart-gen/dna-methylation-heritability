#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=regctx_prox
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/regctx_prox.%A_%a.log
#SBATCH --array=0-11

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

RUNS=(
  "caudate AA shared"
  "caudate EA shared"
  "caudate AA AA_only"
  "caudate EA EA_only"
  "dlpfc AA shared"
  "dlpfc EA shared"
  "dlpfc AA AA_only"
  "dlpfc EA EA_only"
  "hippocampus AA shared"
  "hippocampus EA shared"
  "hippocampus AA AA_only"
  "hippocampus EA EA_only"
)
read -r TISSUE POP VMR_SET <<< "${RUNS[$SLURM_ARRAY_TASK_ID]}"
Rscript "../_h/04.feature_proximity.R" \
  all_individuals "$TISSUE" "$POP" "$VMR_SET"

conda deactivate
log_message "**** Job ends ****"
