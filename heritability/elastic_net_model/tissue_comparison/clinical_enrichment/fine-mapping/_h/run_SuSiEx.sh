#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=1G                # Memory limit
#SBATCH --job-name=run_SuSiEx  # Job name
#SBATCH --output=output_%j.log  # Standard output log
#SBATCH --error=error_%j.log    # Standard error log

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

cd ../SuSiEx/examples
../bin/SuSiEx \
	--sst_file=AFR.sumstats.txt \
	--n_gwas=50000 \
	--ref_file=../../TOPMed_LIBD.AA.998311_998546 \
	--ld_file=../TOPMed_LIBD.AA.998311_998546 \
	--out_dir=../../_m \
	--out_name=TOPMed_LIBD.AA.998311_998546.output \
	--level=0.95 \
	--pval_thresh=1e-5 \
	--maf=0.005 \
	--chr=1 \
	--bp=998311,998546 \
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

log_message "**** Job ends ****"