#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=gengpu
#SBATCH --gres=gpu:a100:2
#SBATCH --job-name=simu_test_500n
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --output=logs/simu_test_500n_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=10:00:00

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

export NUM_SAMPLES=500
export CUPY_CACHE_DIR=/projects/p32505/opt/cupy_cache
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh

log_message "**** Run elastic net ****"
conda activate /projects/p32505/opt/envs/ml

export LD_LIBRARY_PATH=/projects/p32505/opt/envs/ml/lib:$LD_LIBRARY_PATH
SCRATCH_BASE="/scratch/$USER"
JOB_TMP="${SCRATCH_BASE}/tmp_${SLURM_JOB_ID:-interactive}"
mkdir -p "$JOB_TMP"

export TMPDIR="$JOB_TMP"
export DASK_TEMPORARY_DIRECTORY="$TMPDIR"
export DASK_DISTRIBUTED__WORKER__TEMPORARY_DIRECTORY="$TMPDIR"

# Helper on HPC
export UCX_TLS=tcp
export DASK_DISTRIBUTED__COMM__TIMEOUTS__CONNECT="60s"
export DASK_DISTRIBUTED__COMM__TIMEOUTS__TCP="60s"

python ../_h/boosting_elastic_net.py
conda deactivate

log_message "**** Job ends ****"
