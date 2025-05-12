#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=00:10:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=5G                # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=combine_summary # Job name
#SBATCH --output=logs/combine_summary.%j.log # Standard output log

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

## List current modules for reproducibility

module purge
module list

# Set path variables 
OUTPUT="./greml_summary_combined.tsv"
INDEX=1

# Run job
log_message "**** Combining summary results for all brain regions ****"
find ../../../*/_m/summary/greml_summary.tsv -type f | while read -r file; do
    # Extract REGION from path — assumes pattern: gcta/REGION/_m/summary/greml_summary.tsv
    REGION=$(basename "$(dirname "$(dirname "$(dirname "$file")")")")

    if [[ $INDEX -eq 1 ]]; then
        # Write header with REGION column
        head -n 1 "$file" | awk -F'\t' 'BEGIN{OFS="\t"} {print $0, "Region"}' > "$OUTPUT"
        INDEX=0
    fi

    # Append data with REGION column
    tail -n +2 "$file" | awk -F'\t' -v r="$REGION" 'BEGIN{OFS="\t"} {print $0, r}' >> "$OUTPUT"
done


log_message "**** Job ends ****"
