#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=1G                # Memory limit
#SBATCH --job-name=liftOver  # Job name
#SBATCH --output=logs/output_%j.log  # Standard output log
#SBATCH --error=logs/error_%j.log    # Standard error log

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

INPUT_FILE=/projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/fine-mapping/SuSiEx/examples/AFR.sumstats.txt

echo "**** Preparing input for liftOver ****"
awk 'BEGIN{OFS="\t"} NR>1 {print "chr"$1, $3-1, $3, $2, $9, "."}' $INPUT_FILE > AFR.bed

echo "**** Running liftOver ****"
conda activate /projects/p32505/opt/envs/epigenomics
/projects/p32505/opt/envs/epigenomics/lib/R/bin/Rscript ../_h/liftOver.R
conda deactivate

echo "**** Processing liftOver output ****"
awk 'BEGIN{OFS="\t"} {print $4, $3}' AFR_lifted.bed > snp2bp.txt

awk 'BEGIN{FS=OFS="\t"}
NR==FNR {bp[$1]=$2; next} 
NR!=FNR {
    if(FNR==1){print; next}  # keep header
    if($2 in bp) $3=bp[$2];  # replace bp with lifted value
    print
}' snp2bp.txt $INPUT_FILE > AFR_lifted.sumstats.txt

log_message "**** Job ends ****"