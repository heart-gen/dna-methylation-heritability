#!/bin/bash

#  After running 'install_software.sh', this should point to the directory
#  where this repo was cloned, and not say "$PWD"

#SBATCH --account=b1042        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G                # Memory limit
#SBATCH --job-name=BiocMAP
#SBATCH -o ./00_run_first_half_slurm.log
#SBATCH -e ./00_run_first_half_slurm.log

log_message "**** Job starts ****"

echo "**** QUEST info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility

module purge
module load fastqc/0.12.0
module list

# Your commands go here

ORIG_DIR=$PWD

export _JAVA_OPTIONS="-Xms8g -Xmx10g"

$ORIG_DIR/Software/bin/nextflow $ORIG_DIR/first_half.nf \
    --annotation "$ORIG_DIR/ref" \
    --sample "paired" \
    --reference "hg38" \
    --all_alignments \
    --input /projects/p32505/users/elisa/projects/dna-methylation-heritability/testing/dlpfc/_h/BiocMAP/benchmark/samples.manifest \
    --output "./results" \
    --trim_mode "adaptive" \
    -profile first_half_slurm
