#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=10G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=1-9974%250
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=greml     # Job name
#SBATCH --output=/dev/null      # Standard output log
#SBATCH --error=/dev/null       # Standard error log

## Edit with your job command
REGION_LIST="./vmr_list.txt"
CHR_FILE="/projects/b1213/resources/genomes/human/gencode-v47/fasta/chromosome_sizes.txt"
DATA="/projects/p32505/users/alexis/projects/dna-methylation-heritability/inputs/genotypes"
WORKING="./"
OUTPUT="./h2"

ENV_PATH="/projects/p32505/opt/env"

# Get the current sample name from the sample list
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
exec > >(tee -a "$LOG_DIR/greml_${SLURM_ARRAY_TASK_ID}_out.log")
exec 2> >(tee -a "$LOG_DIR/greml_${SLURM_ARRAY_TASK_ID}_err.log" >&2)

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

##### GREML-LDMS #####
# Check if SNP data exists
BFILE="$WORKING/plink_format/chr_${CHR}/TOPMed_LIBD.AA.${START}_${END}"

if [ ! -f "${BFILE}.bed" ] || [ ! -f "${BFILE}.bim" ] || [ ! -f "${BFILE}.fam" ]; then
    log_message "SNP files for region $CHR: $START-$END not found. Skipping."
    exit 0
fi

echo "Performing GREML analysis on $CHR: $START-$END"

# Calculate SNP LD scores
gcta64 --bfile $WORKING/plink_format/chr_${CHR}/TOPMed_LIBD.AA.${START}_${END} \
       --ld-score-region 200 \
       --out $CHR_DIR/TOPMed_LIBD.AA.${START}_${END}

## Activate conda environment
# Stratify SNPs based on individual LD scores
conda run -p $ENV_PATH/r_env Rscript ../_h/05.stratify_LD.R $CHR $START $END

if [ $? -ne 0 ]; then
    log_message "Error: Conda or script execution failed"
    exit 1
fi

# Make GRM for each group
for i in {1..4} ; do
    gcta64 --bfile $WORKING/plink_format/chr_${CHR}/subset_TOPMed_LIBD.AA.${START}_${END} \
           --extract $CHR_DIR/${START}_${END}_snp_group_${i}.txt \
           --make-grm \
           --out $CHR_DIR/TOPMed_LIBD.AA.${START}_${END}_group_${i}

    echo "$CHR_DIR/TOPMed_LIBD.AA.${START}_${END}_group_${i}" >> \
         $CHR_DIR/${START}_${END}_multi_GRMs.txt
done

# GREML with multiple GRM
gcta64 --reml \
       --mgrm $CHR_DIR/${START}_${END}_multi_GRMs.txt \
       --pheno $WORKING/vmr/chr_${CHR}/${START}_${END}_meth.phen \
       --covar $WORKING/covs/chr_${CHR}/TOPMed_LIBD.AA.covar \
       --qcovar $WORKING/covs/chr_${CHR}/TOPMed_LIBD.AA.qcovar \
       --out $CHR_DIR/TOPMed_LIBD.AA.${START}_${END}

log_message "**** Job ends ****"
