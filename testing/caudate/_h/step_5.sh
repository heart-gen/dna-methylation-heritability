#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=10G               # Memory limit
#SBATCH --mail-type=FAIL
#SBATCH --array=1-12078%250
#SBATCH --mail-user=alexis.bennett@northwestern.edu
#SBATCH --job-name=extract_snp  # Job name
#SBATCH --output=logs/extract_snp_%j_out.log  # Standard output log
#SBATCH --error=logs/extract_snp_%j_err.log    # Standard error log

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
REGION_LIST="./vmr_list.txt"
CHR_FILE="/projects/b1213/resources/genomes/human/gencode-v47/fasta/chromosome_sizes.txt"
DATA="/projects/p32505/projects/dna-methylation-heritability/inputs/genotypes"
OUTPUT="/projects/p32505/users/alexis/projects/dna-methylation-heritability/testing/caudate/_m"

ENV="/projects/p32505/opt/env"

# Get the current sample name from the sample list
REGION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $REGION_LIST)
CHR=$(echo "$REGION" | awk '{print $1}')
START=$(echo "$REGION" | awk '{print $2}')
END=$(echo "$REGION" | awk '{print $3}')

echo "Perfoming GREML analysis on $CHR: $START-$END"

##### GREML-LDMS #####
gcta64 --bfile $OUTPUT/chr_$CHR/TOPMed_LIBD.AA.${START}_${END} --ld-score-region 200 --out $OUTPUT/chr_$CHR/TOPMed_LIBD.AA.${START}_${END}

# Activating conda environment
$ENV/R_env/bin/Rscript ../_h/05.stratify_LD.R

# Make GRM for each group
for i in 1:4 ; do
gcta64 --bfile $OUTPUT/chr_$CHR/TOPMed_LIBD.AA.${START}_${END} --extract $OUTPUT/chr_$CHR/${START}_${END}_snp_group_${i}.txt --make-grm --out $OUTPUT/chr_$CHR/TOPMed_LIBD.AA.${START}_${END}_group_${i}

echo "$OUTPUT/chr_$CHR/TOPMed_LIBD.AA.${START}_${END}_group_${i}" >> $OUTPUT/chr_$CHR/${START}_${END}_multi_GRMs.txt
done

# GREML with multiple GRM
gcta64 --reml --mgrm $OUTPUT/chr_$CHR/${START}_${END}_multi_GRMs.txt --pheno $OUTPUT/chr_$CHR/${START}_${END}_meth.phen --covar $OUTPUT/chr_$CHR/TOPMed_LIBD.AA.covar --qcovar $OUTPUT/chr_$CHR/TOPMed_LIBD.AA.qcovar --out TOPMed_LIBD.AA.${START}_${END}
