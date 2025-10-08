#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=16G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=1-9576%250
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=extract_snp  # Job name
#SBATCH --output=/dev/null      # Standard output log
#SBATCH --error=/dev/null       # Standard error log

## Edit with your job command
REGION_LIST="./vmr.bed"
SAMPLE_LIST="./samples.txt"
CHR_FILE="/projects/b1213/resources/genomes/human/gencode-v47/fasta/chromosome_sizes.txt"
DATA="/projects/b1213/users/alexis/projects/dna-methylation-heritability/inputs/genotypes"
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
END_POS=$((END + WINDOW))

echo "Extracting SNPs from all subjects on $CHR: $START-$END ($WINDOW bp window)"

if (( START_POS <= 0 )); then
    echo "ERROR: Start position is below zero."
    exit 1
fi

if (( END_POS >= CHR_SIZE )); then
    echo "ERROR: End position exceeds Chromosome $CHR size."
    exit 1
fi

plink2 --pfile "$DATA/TOPMed_LIBD.AA" \
       --chr "$CHR" \
       --from-bp "$START_POS" \
       --to-bp "$END_POS" \
       --make-bed \
       --no-parents \
       --no-sex \
       --no-pheno \
       --out "$CHR_DIR/TOPMed_LIBD.AA.${START}_${END}"

echo "Extracting SNPs from AA subjects on $CHR: $START-$END ($WINDOW bp window)"

# Subset of SNPs in AA cohort
plink2 --pfile "$DATA/TOPMed_LIBD.AA" \
       --chr "$CHR" \
       --from-bp "$START_POS" \
       --to-bp "$END_POS" \
       --make-bed \
       --keep "$SAMPLE_LIST" \
       --no-parents \
       --no-sex \
       --no-pheno \
       --out "$CHR_DIR/subset_TOPMed_LIBD.AA.${START}_${END}"

log_message "**** Job ends ****"
