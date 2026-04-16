#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=nh_subgroup_classify
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --output=logs/subgroup_classify_%A.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=00:30:00

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

module purge
module list

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh

mkdir -p logs ../_m

log_message "**** Step 1: Classify AA non-heritable VMRs by EA heritability status ****"
conda activate /projects/p32505/opt/envs/genomics
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
python 01.subgroup_classification.py
conda deactivate

log_message "**** Step 2: Fisher's exact enrichment per subgroup ****"
conda activate /projects/p32505/opt/envs/genomics
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
python 02.fishers_enrichment.py
conda deactivate

log_message "**** Job ends ****"
