#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_cpg_phenotypes
#SBATCH --output=logs/cpg_phenotypes.%A_%a.log

# CPU: build per-chromosome CpG BEDs within VMRs, then merge autosomal BED
# Requires step_1 outputs. Submit from _m/:
#   sbatch ../_h/step_2.sh

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
echo "**** QUEST info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Node name: ${SLURM_NODENAME:-N/A}"
echo "Array task: ${SLURM_ARRAY_TASK_ID:-N/A}"

module purge
module list
mkdir -p logs

# bgzip/tabix for 02b merge (not present on bare node PATH)
log_message "**** Loading conda environment ****"
# shellcheck disable=SC1091
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics
command -v bgzip >/dev/null && command -v tabix >/dev/null \
  || { log_message "ERROR: bgzip/tabix missing after activating genomics env"; exit 1; }

REGIONS=(caudate dlpfc hippocampus)
if [[ -z "${REGION:-}" ]]; then
  REGION="${REGIONS[${SLURM_ARRAY_TASK_ID}]}"
fi

PREP="../${REGION}/_m/prepared"
CPG_ROOT="/projects/b1213/users/alexis/projects/dna-methylation-heritability/vmr-analysis/${REGION}/_m/cpg"

log_message "Building CpG phenotype BEDs for ${REGION}"
for d in "${CPG_ROOT}"/chr_*; do
  chrom="${d##*_}"
  if [[ "${chrom}" =~ ^(X|Y|M|MT)$ ]]; then
    continue
  fi
  log_message "  chr${chrom}"
  python3 "../_h/02a_prepare_cpg_bed.py" --region "${REGION}" --chrom "${chrom}"
done

log_message "Merging autosomal CpG BEDs for ${REGION}"
python3 "../_h/02b_merge_cpg_beds.py" --prepared-dir "${PREP}"

conda deactivate
log_message "**** Job ends ****"
