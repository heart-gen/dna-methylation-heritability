#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=simu_array
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/simulation_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=150G
#SBATCH --array=0-2

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

echo "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME:-N/A}"
echo "Hostname: ${HOSTNAME}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-N/A}"

module purge
module list

log_message "**** Run simulation ****"
ENV_PATH="/projects/p32505/opt/envs"

# Sample sizes corresponding to array index
LD_DECAY=(0.5 0.6 0.7)
CURRENT_DECAY=${LD_DECAY[$SLURM_ARRAY_TASK_ID]}

conda run -p "${ENV_PATH}/genomics" \
      python ../_h/01.simulated-data.py \
      --num_phenotypes 1000 \
      --num_samples 200 \
      --ld_decay "${CURRENT_DECAY}" \
      --output_dir ./ld_${CURRENT_DECAY}_sim_200_indiv

log_message "**** Job ends ****"
