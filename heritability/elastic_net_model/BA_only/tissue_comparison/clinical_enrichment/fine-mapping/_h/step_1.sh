#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=20gb
#SBATCH --job-name=extract_rsids
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --output=logs/extract_rsids.%j.log
#SBATCH --array=1-22

# Function to echo with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

VMR_DIR="/projects/b1213/users/alexis/projects/dna-methylation-heritability/vmr-analysis/hippocampus/_m/plink_format/chr_$SLURM_ARRAY_TASK_ID"
OUT_DIR="./snplists/hippocampus/chr_${SLURM_ARRAY_TASK_ID}"

mkdir -p ${OUT_DIR}

for BIM in ${VMR_DIR}/*.bim; do

    PREFIX=$(basename ${BIM} .bim)

    echo "Extracting rsIDs for ${PREFIX}"

    awk '{print $2}' ${BIM} | awk -F'_' '{print $NF}' \
        > ${OUT_DIR}/${PREFIX}_rsids.txt

done

echo "All rsID extraction complete."

log_message "**** Job ends ****"