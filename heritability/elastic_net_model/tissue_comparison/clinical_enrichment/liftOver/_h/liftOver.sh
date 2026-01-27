#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=5G                # Memory limit
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

INPUT_FILE=/projects/b1213/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/fine-mapping/SuSiEx/examples/AFR.sumstats.txt

awk 'BEGIN{OFS="\t"} NR>1 {print "chr"$1, $3-1, $3, $2}' $INPUT_FILE \
| sort -k1,1 -k2,2n > gwas.bed

# Runs liftOver
conda activate /projects/p32505/opt/envs/epigenomics
/projects/p32505/opt/envs/epigenomics/lib/R/bin/Rscript ../_h/liftOver.R
conda deactivate

awk '{print $4, $3}' gwas_lifted.bed > lifted.map

awk '
NR==FNR {
    newpos[$1]=$2
    next
}
NR==1 {
    print
    next
}
{
    split($2, a, ":")     # a[1]=id, a[2]=old bp, a[3]=A1, a[4]=A2
    bp38 = newpos[$2]

    if (bp38 != "") {
        $2 = a[1] ":" bp38 ":" a[3] ":" a[4]  # new SNP ID
        $3 = bp38                              # new bp column
    }
    print
}
' OFS="\t" lifted.map $INPUT_FILE > gwas.sumstats.txt

rm gwas.bed gwas_lifted.bed lifted.map

log_message "**** Job ends ****"