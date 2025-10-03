#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=gengpu
#SBATCH --gres=gpu:a100:1
#SBATCH --job-name=boosting_elastic_h2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --output=logs/boosting_elastic_h2_%A_%a.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=01:00:00

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
echo "OFFSET: ${OFFSET}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"

module purge
module load cuda/12.0.1-gcc-12.3.0
module list

export NUM_SAMPLES=100
export CUPY_CACHE_DIR=/projects/p32505/opt/cupy_cache
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh

log_message "**** Run elastic net ****"
conda activate /projects/p32505/opt/envs/ml
export LD_LIBRARY_PATH=/projects/p32505/opt/envs/ml/lib:$LD_LIBRARY_PATH
python ../_h/boosting_elastic_net.py
conda deactivate
#ENV_PATH="/projects/p32505/opt/env"
#conda run -p "${ENV_PATH}/AI_env" python ../_h/boosting_elastic_net.py

log_message "**** Job ends ****"
