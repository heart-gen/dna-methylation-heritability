#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=16G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=0-5999%250
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
SAMPLE_SIZES=(100 150 200 250 500 1000)
NUM_PHEN=1000

IDX=$(( SLURM_ARRAY_TASK_ID / NUM_PHEN ))
PHEN_INDEX=$((SLURM_ARRAY_TASK_ID % NUM_PHEN))
SAMPLE_SIZE=${SAMPLE_SIZES[$IDX]}

OUTPUT="./h2/sim_${SAMPLE_SIZE}_indiv"
mkdir -p "$OUTPUT"

ENV_PATH="/projects/p32505/opt/env"

##### GREML-LDMS #####
# Check if SNP and LD score data exists
## This checks only for PLINK files, needs to be updated for LD files
BFILE="$SIM/sim_${SAMPLE_SIZE}_indiv/plink_sim/"
required_files=(simulated.bed simulated.bim simulated.fam)
for f in "${required_files[@]}"; do
    if [[ ! -f "${PLINK_DIR}/$f" ]]; then
        log_message "SNP files not found. Skipping."
        exit 0
    fi
done

log_message "Performing GREML analysis on simulated phenotype ${PHEN_INDEX} with ${SAMPLE_SIZE} samples"

# GREML with multiple GRM
gcta64 --reml \
       --mgrm $OUTPUT/sim_${SAMPLE_SIZE}_indiv_multi_GRMs.txt \
       --pheno $SIM/sim_${SAMPLE_SIZE}_indiv/simulated.phen \
       --mpheno $PHEN_INDEX \
       --out $OUTPUT/greml_pheno$PHEN_INDEX

log_message "**** Job ends ****"
