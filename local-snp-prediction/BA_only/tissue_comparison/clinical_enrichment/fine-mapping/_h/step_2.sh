#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=08:00:00
#SBATCH --mem=40gb
#SBATCH --job-name=match_and_subset
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --output=logs/match_and_subset.%j.log
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

GWAS="/projects/b1213/resources/gwas/PD/data/GCST009325.h.tsv.gz"
SNPLISTS="./snplists/hippocampus/chr_$SLURM_ARRAY_TASK_ID"          # where *_rsids.txt currently are
OUT_DIR="./gwas/pd/hippocampus/chr_$SLURM_ARRAY_TASK_ID"     # NEW output directory

mkdir -p ${OUT_DIR}

for RSID_FILE in ${SNPLISTS}/*_rsids.txt; do

    PREFIX=$(basename ${RSID_FILE} _rsids.txt)

    echo "Matching SNPs for ${PREFIX}"

    OUT_GWAS=${OUT_DIR}/${PREFIX}_gwas_subset.tsv
    OUT_SNPS=${OUT_DIR}/${PREFIX}_matched_snps.txt

    zcat ${GWAS} | awk \
        -v snps=${RSID_FILE} \
        -v out_gwas=${OUT_GWAS} \
        -v out_snps=${OUT_SNPS} '
        
        BEGIN {
            # Load SNP list and preserve order
            while ((getline line < snps) > 0) {
                order[line] = ++i
                snplist[i] = line
            }
        }

        NR==1 {
            print $0 > out_gwas   # write header
            next
        }

        {
            rsid = $9
            if (rsid in order) {
                gwas[rsid] = $0
            }
        }

        END {
            for (j = 1; j <= i; j++) {
                if (snplist[j] in gwas) {
                    print gwas[snplist[j]] >> out_gwas
                    print snplist[j] >> out_snps
                }
            }
        }
    '

done

echo "All GWAS matching complete."

log_message "**** Job ends ****"