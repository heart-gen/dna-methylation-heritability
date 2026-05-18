#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=80G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=regctx_assoc
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/regctx_assoc.%A_%a.log
#SBATCH --array=0-17

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

# Shared matched VMRs: one association run per tissue/modality (population AA).
# EA/shared is skipped in 01.run_local_associations.R; steps 02–03 reuse AA inputs.
RUNS=(
  "caudate expression AA 250000 shared"
  "caudate psi AA 250000 shared"
  "caudate expression AA 250000 AA_only"
  "caudate psi AA 250000 AA_only"
  "caudate expression EA 250000 EA_only"
  "caudate psi EA 250000 EA_only"
  "dlpfc expression AA 250000 shared"
  "dlpfc psi AA 250000 shared"
  "dlpfc expression AA 250000 AA_only"
  "dlpfc psi AA 250000 AA_only"
  "dlpfc expression EA 250000 EA_only"
  "dlpfc psi EA 250000 EA_only"
  "hippocampus expression AA 250000 shared"
  "hippocampus psi AA 250000 shared"
  "hippocampus expression AA 250000 AA_only"
  "hippocampus psi AA 250000 AA_only"
  "hippocampus expression EA 250000 EA_only"
  "hippocampus psi EA 250000 EA_only"
)

log_message "Running local associations"

read -r TISSUE MODALITY POP WINDOW VMR_SET <<< "${RUNS[$SLURM_ARRAY_TASK_ID]}"
Rscript "../_h/01.run_local_associations.R" \
  all_individuals "$TISSUE" "$MODALITY" "$POP" "$WINDOW" "$VMR_SET"

conda deactivate
log_message "**** Job ends ****"
