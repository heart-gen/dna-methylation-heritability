#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=16G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=1-10216%250
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=extract_snp  # Job name
#SBATCH --output=/dev/null      # Standard output log
#SBATCH --error=/dev/null       # Standard error log

## Define paths
REGION_LIST="./vmr.bed"
SAMPLE_LIST="./samples.txt"
CHR_FILE="/projects/b1213/resources/genomes/human/gencode-v47/fasta/chromosome_sizes.txt"

# Repository root, resolved the way here::here() does in the R scripts
REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
while [ "$REPO_DIR" != "/" ] && [ ! -d "$REPO_DIR/.git" ]; do
    REPO_DIR=$(dirname "$REPO_DIR")
done

PLINK2="/projects/p32505/opt/bin/plink2"

DATA="$REPO_DIR/inputs/genotypes"
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

## Check inputs and software
##
## plink2 comes from the shared opt tree rather than a module: the genomics
## conda env ships plink 1.9 only, which cannot read .pgen/.pvar input.

if [ ! -d "$REPO_DIR/.git" ]; then
    echo "ERROR: could not locate repository root from ${SLURM_SUBMIT_DIR:-$PWD}"
    exit 1
fi

if [ ! -x "$PLINK2" ]; then
    echo "ERROR: plink2 not found at $PLINK2"
    exit 1
fi

"$PLINK2" --version

# check chromosome size information

WINDOW=500000
CHR_SIZE=$(grep "^chr1[[:space:]]" $CHR_FILE | cut -f2)

START_POS=$((START - WINDOW))
END_POS=$((END + WINDOW))

if (( START_POS <= 0 )); then
    echo "ERROR: Start position is below zero."
    exit 1
fi

if (( END_POS >= CHR_SIZE )); then
    echo "ERROR: End position exceeds Chromosome $CHR size."
    exit 1
fi

echo "Extracting SNPs from AA subjects on $CHR: $START-$END ($WINDOW bp window)" 

# Subset of SNPs in AA cohort
"$PLINK2" --pfile "$DATA/TOPMed_LIBD.AA" \
          --chr "$CHR" \
          --from-bp "$START_POS" \
          --to-bp "$END_POS" \
          --make-bed \
          --keep "$SAMPLE_LIST" \
          --no-parents \
          --no-sex \
          --no-pheno \
          --out "$CHR_DIR/TOPMed_LIBD.AA.${START}_${END}"

log_message "**** Job ends ****"
