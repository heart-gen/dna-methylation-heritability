#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=16G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=1-12078%250
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=extract_snp  # Job name
#SBATCH --output=/dev/null      # Standard output log
#SBATCH --error=/dev/null       # Standard error log

## Edit with your job command
REGION_LIST="./vmr_list.txt"
SAMPLE_LIST="./samples.txt"
CHR_FILE="/projects/b1213/resources/genomes/human/gencode-v47/fasta/chromosome_sizes.txt"
DATA="/projects/p32505/users/alexis/projects/dna-methylation-heritability/inputs/genotypes"
OUTPUT="./plink_format"

# Get the current region name from the region list
REGION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $REGION_LIST)
CHR=$(echo "$REGION" | awk '{print $1}')
START=$(echo "$REGION" | awk '{print $2}')
END=$(echo "$REGION" | awk '{print $3}')

# Create directories for each chr
CHR_DIR="$OUTPUT/chr_${CHR}"
mkdir -p "$CHR_DIR"
LOG_DIR="logs/chr_${CHR}"
mkdir -p "$LOG_DIR"

# Redirect output and error logs to chr-specific log files
exec > >(tee -a "$LOG_DIR/extract_snp_${SLURM_ARRAY_TASK_ID}_out.log")
exec 2> >(tee -a "$LOG_DIR/extract_snp_${SLURM_ARRAY_TASK_ID}_err.log" >&2)

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
module load plink/2.0-alpha-3.3
module list

# check chromosome size information
WINDOW=500000
CHR_SIZE=$(grep "^chr1[[:space:]]" $CHR_FILE | cut -f2)

START_POS=$((START - WINDOW))

if (( START_POS <= 0 )); then
    echo "ERROR: Start position is below zero."
    exit 1
fi

END_POS=$((END + WINDOW))

if (( END_POS >= CHR_SIZE )); then
    echo "ERROR: End position exceeds Chromosome $CHR size."
    exit 1
fi

echo "Extracting SNPs from all subjects on $CHR: $START_POS-$END_POS"

plink2 --pfile "$DATA/TOPMed_LIBD.AA" \
       --chr "$CHR" \
       --from-bp "$START_POS" \
       --to-bp "$END_POS" \
       --make-bed \
       --out "$CHR_DIR/TOPMed_LIBD.AA.${START}_${END}"

echo "Extracting SNPs from AA subjects on $CHR: $START_POS-$END_POS"

# Subset of SNPs in AA cohort
plink2 --pfile "$DATA/TOPMed_LIBD.AA" \
       --chr "$CHR" \
       --from-bp "$START_POS" \
       --to-bp "$END_POS" \
       --keep "$SAMPLE_LIST" \
       --make-bed \
       --out "$CHR_DIR/subset_TOPMed_LIBD.AA.${START}_${END}"

log_message "**** Job ends ****"