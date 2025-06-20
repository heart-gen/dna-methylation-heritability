#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=16G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=0-131
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=calculate_LD     # Job name
#SBATCH --output=logs/calculate_LD_%j.out    # Standard output log

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
module load gcta/1.94.0
module list

## Edit with your job command
SIM="../../../inputs/simulated-data/_m" # Switched to relative PATH for flexibility
ENV_PATH="/projects/p32505/opt/env"

SAMPLE_SIZES=(100 150 200 250 500 1000)
NUM_CHR=22

IDX=$(( SLURM_ARRAY_TASK_ID / NUM_CHR ))
CHR=$(( (SLURM_ARRAY_TASK_ID % NUM_CHR) + 1 ))
SAMPLE_SIZE=${SAMPLE_SIZES[$IDX]}

OUTPUT="./h2/sim_${SAMPLE_SIZE}_indiv"
mkdir -p "$OUTPUT"

##### GREML-LDMS #####
# Check if SNP data exists
BFILE="$SIM/sim_${SAMPLE_SIZE}_indiv/plink_sim/"

if [ ! -f "${BFILE}.bed" ] || [ ! -f "${BFILE}.bim" ] || [ ! -f "${BFILE}.fam" ]; then
    log_message "SNP files not found. Skipping."
    exit 0
fi

log_message "Calaculating SNP LD scores for simulated data with $SAMPLE_SIZE samples"

# Calculate SNP LD scores
gcta64 --bfile $BFILE \
    --chr $CHR \
    --ld-score-region 200 \
    --out $OUTPUT/sim_${SAMPLE_SIZE}_indiv_chr${CHR}

log_message "**** Job ends ****"
