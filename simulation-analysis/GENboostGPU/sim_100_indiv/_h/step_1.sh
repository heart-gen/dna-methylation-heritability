#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=gengpu
#SBATCH --gres=gpu:a100:2
#SBATCH --job-name=simu_test_100n
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --output=logs/simu_test_100n_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=24:00:00

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
module load cuda/12.4.1-gcc-12.3.0
module list

export NUM_SAMPLES=100
export CUPY_CACHE_DIR=/projects/p32505/opt/cupy_cache
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh

log_message "**** Run elastic net ****"
conda activate /projects/p32505/opt/envs/ml
export LD_LIBRARY_PATH=/projects/p32505/opt/envs/ml/lib:$LD_LIBRARY_PATH
python ../_h/simu_test_100n.py
conda deactivate
#ENV_PATH="/projects/p32505/opt/env"
#conda run -p "${ENV_PATH}/AI_env" python ../_h/boosting_elastic_net.py

log_message "**** Job ends ****"