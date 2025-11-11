#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name=region_heritability
#SBATCH --output=logs/00.output_%j.log
#SBATCH --error=logs/00.error_%j.log

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
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

BRAIN_REGIONS=("caudate" "dlpfc" "hippocampus")

for REGION in "${BRAIN_REGIONS[@]}"; do

    PYTHON_SCRIPT="../_h/00.region_heritability.py"
    INPUT_FILE="/projects/b1213/users/alexis/projects/dna-methylation-heritability/heritability/elastic_net_model/${REGION}/_m/${REGION}_summary_elastic-net.tsv"
    OUTPUT_DIR="./vmr/${REGION}"
    mkdir -p "$OUTPUT_DIR"

    echo "Running separation script for region: $REGION"
    python "$PYTHON_SCRIPT" \
        --input_file "$INPUT_FILE" \
        --output_dir "$OUTPUT_DIR"

    echo "Separation completed. Check files in $OUTPUT_DIR"

done
log_message "**** Job ends ****"