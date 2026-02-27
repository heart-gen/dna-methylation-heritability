#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --mem=20gb
#SBATCH --job-name=regenerate_plink
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --output=logs/regenerate_plink.%j.log

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

PLINK_DIR="/projects/b1213/users/alexis/projects/dna-methylation-heritability/vmr-analysis/hippocampus/_m/plink_format/chr_21"
SUBSET_DIR="./pd/chr_21"
OUT_DIR="./ld_matrices/hippocampus/chr_21"

for MATCHED in ${SUBSET_DIR}/*_matched_snps.txt; do

    PREFIX=$(basename ${MATCHED} _matched_snps.txt)
    BIM=${PLINK_DIR}/${PREFIX}.bim
    FULL_PREFIX=${PLINK_DIR}/${PREFIX}

    echo "Regenerating PLINK files for ${PREFIX}"

    EXTRACT_FILE=${SUBSET_DIR}/${PREFIX}_extract.txt

    # Map rsID back to original PLINK SNP ID
    awk -v snps=${MATCHED} '
        BEGIN {
            while ((getline line < snps) > 0) {
                keep[line]=1
            }
        }
        {
            split($2, arr, "_")
            rsid=arr[length(arr)]
            if (rsid in keep) {
                print $2
            }
        }
    ' ${BIM} > ${EXTRACT_FILE}

    # Regenerate PLINK dataset
    ${PLINK} \
        --bfile ${FULL_PREFIX} \
        --extract ${EXTRACT_FILE} \
        --make-bed \
        --out ${OUT_DIR}/${PREFIX}_final

done

echo "All PLINK regeneration complete."

log_message "**** Job ends ****"