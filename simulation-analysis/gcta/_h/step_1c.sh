#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=02:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=16G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=0-7999%250
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=greml     # Job name
#SBATCH --output=logs/greml_%j.out    # Standard output log

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
SIM="../../../inputs/simulated-data/_m"
SAMPLE_SIZES=(100 150 200 250 500 1000 5000 10000)
NUM_PHEN=1000

IDX=$(( SLURM_ARRAY_TASK_ID / NUM_PHEN ))
PHEN_INDEX=$(( (SLURM_ARRAY_TASK_ID % NUM_PHEN) + 1 ))
SAMPLE_SIZE=${SAMPLE_SIZES[$IDX]}

OUTPUT="./h2/sim_${SAMPLE_SIZE}_indiv"
mkdir -p "$OUTPUT"

##### GREML-LDMS #####
# Check if SNP and LD score data exists
PLINK_DIR="$SIM/sim_${SAMPLE_SIZE}_indiv/plink_sim/"
required_files=(simulated.bed simulated.bim simulated.fam)
for f in "${required_files[@]}"; do
    if [[ ! -f "${PLINK_DIR}/$f" ]]; then
        log_message "SNP files not found. Skipping."
        exit 0
    fi
done

LDSCORE="$OUTPUT/sim_${SAMPLE_SIZE}_indiv.score.ld"
if [[ ! -f "$LDSCORE" ]]; then
    log_message "LD score file not found. Skipping"
    exit 0
fi

log_message "Performing GREML analysis on simulated phenotype ${PHEN_INDEX} with ${SAMPLE_SIZE} samples"

# GREML with multiple GRM
gcta64 --reml \
       --mgrm $OUTPUT/sim_${SAMPLE_SIZE}_indiv_multi_GRMs.txt \
       --pheno $SIM/sim_${SAMPLE_SIZE}_indiv/simulated.phen \
       --mpheno $PHEN_INDEX \
       --out $OUTPUT/greml_pheno_$PHEN_INDEX

log_message "**** Job ends ****"
