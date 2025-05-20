#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=1-22
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --job-name=BiocMAP  # Job name
#SBATCH -q shared
#SBATCH -o ./run_first_half_jhpce.log
#SBATCH -e ./run_first_half_jhpce.log

# Log function
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
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## List current modules for reproducibility

module purge
module list

#  After running 'install_software.sh', this should point to the directory
#  where BiocMAP was installed, and not say "$PWD"
ORIG_DIR=/home/uyb5533

module load nextflow/20.01.0
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

nextflow $ORIG_DIR/first_half.nf \
    --sample "paired" \
    --reference "hg38" \
    --input "/users/neagles/wgbs_test" \
    --output "/users/neagles/wgbs_test/out" \
    -w "/fastscratch/myscratch/neagles/nextflow_work" \
    --trim_mode "force" \
    -profile first_half_jhpce
