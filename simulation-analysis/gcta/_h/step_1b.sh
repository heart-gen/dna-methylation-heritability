#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=02:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=16G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --array=0-7
#SBATCH --job-name=stratify_LD     # Job name
#SBATCH --output=logs/stratify_LD_%j.out    # Standard output log

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
SAMPLE_SIZE=${SAMPLE_SIZES[$SLURM_ARRAY_TASK_ID]}

OUTPUT="./h2/sim_${SAMPLE_SIZE}_indiv"
mkdir -p "$OUTPUT"

ENV_PATH="/projects/p32505/opt/env"

# Check for LD score file for each chr
MISSING=0
for CHR in {1..22}; do
    LDSCORE="$OUTPUT/sim_${SAMPLE_SIZE}_indiv_chr${CHR}.score.ld"
    if [[ ! -f "$LDSCORE" ]]; then
        log_message "Missing LD score file: $LDSCORE"
        ((MISSING++))
    fi
done

if [[ $MISSING -ne 0 ]]; then
    log_message "Error: $MISSING chromosome LD score files are missing. Cannot combine."
    exit 1
fi

# Combine chr LD scores to one file
log_message "Combining chromosome LD scores"
awk 'FNR==1 && NR!=1 { next } { print }' $OUTPUT/sim_${SAMPLE_SIZE}_indiv_chr*.score.ld > $OUTPUT/sim_${SAMPLE_SIZE}_indiv.score.ld

export TMPDIR="/projects/p32505/users/alexis/tmp"
mkdir -p "$TMPDIR"

## Activate conda environment
# Stratify SNPs based on individual LD scores
conda run -p $ENV_PATH/r_env Rscript ../_h/01.stratify_LD.R $SAMPLE_SIZE

if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
fi

# Make GRM for each group
for i in {1..4} ; do
    gcta64 --bfile $SIM/sim_${SAMPLE_SIZE}_indiv/plink_sim/simulated \
           --extract $OUTPUT/snp_group_${i}.txt \
           --make-grm \
           --out $OUTPUT/sim_${SAMPLE_SIZE}_indiv_group_${i}

    echo "$OUTPUT/sim_${SAMPLE_SIZE}_indiv_group_${i}" >> \
         $OUTPUT/sim_${SAMPLE_SIZE}_indiv_multi_GRMs.txt
done

log_message "**** Job ends ****"
