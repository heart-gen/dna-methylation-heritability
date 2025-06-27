## Load the package at the top of your script
library("sessioninfo")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

#!/bin/bash

#  After running 'install_software.sh', this should point to the directory
#  where this repo was cloned, and not say "$PWD"

#SBATCH --account=b1042       
#SBATCH --partition=short       
#SBATCH --time=01:00:00         
#SBATCH --nodes=1               
#SBATCH --ntasks-per-node=1     
#SBATCH --mem=25G                
#SBATCH --job-name=BiocMAP
#SBATCH -o ./01_run_second_half_slurm.log
#SBATCH -e ./01_run_second_half_slurm.log

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
module load nextflow/22.04.4
module list
module load java/jdk11.0.10

# Your commands go here

ORIG_DIR=$PWD

export _JAVA_OPTIONS="-Xms8g -Xmx10g"

$ORIG_DIR/Software/bin/nextflow $ORIG_DIR/second_half.nf \
    --annotation "$ORIG_DIR/ref" \
    --sample "paired" \
    --reference "hg38" \
    --input /projects/p32505/users/elisa/projects/dna-methylation-heritability/testing/dlpfc/_h/BiocMAP/benchmark/rules.txt \
    --output "./results" \
    --use_bme \
    -profile second_half_slurm
