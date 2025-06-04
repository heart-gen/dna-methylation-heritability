#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=normal
#SBATCH --job-name=simu_array
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/simulation_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=150G
#SBATCH --array=0-5

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

echo "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"

module purge
module list

log_message "**** Run simulation ****"
ENV_PATH="/projects/p32505/opt/env"

# Sample sizes corresponding to array index
SAMPLES=(100 150 200 250 500 1000)
CURRENT_SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

conda run -p "${ENV_PATH}/AI_env" \
      python ../_h/01.simulated-data.py \
      --num_phenotypes 1000 \
      --num_samples $CURRENT_SAMPLE \
      --ld_decay 0.80 \
      --output_dir ./sim_${CURRENT_SAMPLE}_indiv

log_message "**** Job ends ****"
