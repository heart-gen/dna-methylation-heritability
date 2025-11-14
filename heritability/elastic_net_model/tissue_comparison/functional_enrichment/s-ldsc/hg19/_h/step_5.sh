#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=combine_files
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --output=logs/05.output_combine_files.log
#SBATCH --error=logs/05.error_combine_files.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=00:10:00

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

#!/bin/bash

# Define paths to all LDSC output files
# Format: disorder_region_h2.txt
LDSC_FILES=(
"results/ad/caudate/heritable_hg19/ad_caudate_heritable_hg19.results"
"results/ad/caudate/non_heritable_hg19/ad_caudate_non_heritable_hg19.results"
"results/ad/dlpfc/heritable_hg19/ad_dlpfc_heritable_hg19.results"
"results/ad/dlpfc/non_heritable_hg19/ad_dlpfc_non_heritable_hg19.results"
"results/ad/hippocampus/heritable_hg19/ad_hippocampus_heritable_hg19.results"
"results/ad/hippocampus/non_heritable_hg19/ad_hippocampus_non_heritable_hg19.results"
"results/mdd/caudate/heritable_hg19/mdd_caudate_heritable_hg19.results"
"results/mdd/caudate/non_heritable_hg19/mdd_caudate_non_heritable_hg19.results"
"results/mdd/dlpfc/heritable_hg19/mdd_dlpfc_heritable_hg19.results"
"results/mdd/dlpfc/non_heritable_hg19/mdd_dlpfc_non_heritable_hg19.results"
"results/mdd/hippocampus/heritable_hg19/mdd_hippocampus_heritable_hg19.results"
"results/mdd/hippocampus/non_heritable_hg19/mdd_hippocampus_non_heritable_hg19.results"
"results/scz/caudate/heritable_hg19/scz_caudate_heritable_hg19.results"
"results/scz/caudate/non_heritable_hg19/scz_caudate_non_heritable_hg19.results"
"results/scz/dlpfc/heritable_hg19/scz_dlpfc_heritable_hg19.results"
"results/scz/dlpfc/non_heritable_hg19/scz_dlpfc_non_heritable_hg19.results"
"results/scz/hippocampus/heritable_hg19/scz_hippocampus_heritable_hg19.results"
"results/scz/hippocampus/non_heritable_hg19/scz_hippocampus_non_heritable_hg19.results"
# Repeat for the other h2 metric files
)

# Join all paths with a comma (for Python to read)
echo "${LDSC_FILES[@]}" > ldsc_file_list.txt

python ../_h/05.make_plots.py

log_message "**** Job ends ****"