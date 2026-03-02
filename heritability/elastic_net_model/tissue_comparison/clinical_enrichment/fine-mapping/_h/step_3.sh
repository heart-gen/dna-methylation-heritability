#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --mem=20gb
#SBATCH --job-name=regenerate_plink
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --output=logs/regenerate_plink.%j.log
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

module load plink/1.9

PLINK_DIR="/projects/b1213/users/alexis/projects/dna-methylation-heritability/vmr-analysis/hippocampus/_m/plink_format/chr_$SLURM_ARRAY_TASK_ID"
SUBSET_DIR="./gwas/pd/hippocampus/chr_$SLURM_ARRAY_TASK_ID"
OUT_DIR="./plink/hippocampus/chr_$SLURM_ARRAY_TASK_ID"
mkdir -p ${OUT_DIR}

for MATCHED in ${SUBSET_DIR}/*_matched_snps.txt; do

    PREFIX=$(basename ${MATCHED} _matched_snps.txt)
    FULL_PREFIX=${PLINK_DIR}/${PREFIX}
    BIM=${PLINK_DIR}/${PREFIX}.bim

    echo "Regenerating PLINK files for ${PREFIX} (rsID-only IDs)"

    EXTRACT_FILE=${OUT_DIR}/${PREFIX}_extract_original_ids.txt
    RENAME_FILE=${OUT_DIR}/${PREFIX}_rename.txt

    # Step 1: Map rsID -> original SNP ID (for extraction)
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

    # Step 2: Create rename file (original_ID  rsID)
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
                print $2, rsid
            }
        }
    ' ${BIM} > ${RENAME_FILE}

    # Step 3: Extract matched SNPs
    plink \
        --bfile ${FULL_PREFIX} \
        --extract ${EXTRACT_FILE} \
        --make-bed \
        --out ${OUT_DIR}/${PREFIX}_tmp

    # Step 4: Rename SNP IDs to rsIDs only
    plink \
        --bfile ${OUT_DIR}/${PREFIX}_tmp \
        --update-name ${RENAME_FILE} \
        --make-bed \
        --out ${OUT_DIR}/${PREFIX}_final

    # Optional: remove temporary files
    rm ${OUT_DIR}/${PREFIX}_tmp.*

    echo "Finished ${PREFIX}"

done

echo "All PLINK regeneration complete."

log_message "**** Job ends ****"