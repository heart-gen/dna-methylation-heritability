#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=02:30:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G                # Memory limit
#SBATCH --job-name=run_SuSiEx  # Job name
#SBATCH --output=logs/output_%j.log  # Standard output log
#SBATCH --error=logs/error_%j.log    # Standard error log
#SBATCH --array=1-22          # Array job range

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

# Create output directory
OUT_DIR=results/hippocampus/chr_$SLURM_ARRAY_TASK_ID
mkdir -p $OUT_DIR

# Define chromosome directory and reference directory
CHR_DIR="./ld_matrices/hippocampus/chr_${SLURM_ARRAY_TASK_ID}"

for dir in "$CHR_DIR"/*/; do
    dir_name=$(basename "$dir")

    region=($(echo "$dir_name" | grep -oE '[0-9]+'))

    if [ "${#region[@]}" -eq 2 ]; then
        START=${region[0]}
        END=${region[1]}
        echo "$dir_name: START=$START, END=$END"
    else
        echo "$dir_name: numbers found = ${region[*]}"
    fi

	START_POS=$((START - 500000))
	END_POS=$((END + 500000))

	../SuSiEx/bin/SuSiEx \
		--sst_file=../MDD_harmonized.sumstats.txt \
		--n_gwas=50000 \
		--ref_file=./ld_matrices/hippocampus/chr_${SLURM_ARRAY_TASK_ID}/${START}_${END}/TOPMed_LIBD.AA.${START}_${END} \
		--ld_file=./ld_matrices/hippocampus/chr_${SLURM_ARRAY_TASK_ID}/${START}_${END}/TOPMed_LIBD.AA.${START}_${END} \
		--out_dir=$OUT_DIR \
		--out_name=TOPMed_LIBD.AA.${START}_${END}.output \
		--level=0.95 \
		--pval_thresh=1e-5 \
		--maf=0.005 \
		--chr=${SLURM_ARRAY_TASK_ID} \
		--bp=$START_POS,$END_POS \
		--snp_col=2 \
		--chr_col=1 \
		--bp_col=3 \
		--a1_col=4 \
		--a2_col=5 \
		--eff_col=6 \
		--se_col=7 \
		--pval_col=9 \
		--plink=../utilities/plink \
		--mult-step=True \
		--keep-ambig=True \
		--threads=16
done

log_message "**** Job ends ****"