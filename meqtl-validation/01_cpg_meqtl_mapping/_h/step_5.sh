#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --array=0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_cpg_qc
#SBATCH --output=logs/cpg_qc.%A_%a.log

# CPU: QC summaries from TensorQTL cis results
# Requires step_4. Submit from _m/:
#   sbatch ../_h/step_5.sh

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

log_message "**** Loading conda environment ****"
# shellcheck disable=SC1091
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

REGIONS=(caudate dlpfc hippocampus)
if [[ -z "${REGION:-}" ]]; then
  REGION="${REGIONS[${SLURM_ARRAY_TASK_ID}]}"
fi

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
# Submit from _m/; do not resolve $0 (SLURM copies the batch script to spool).
H="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h"
CIS="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m/tensorqtl/cpg_meqtl_${REGION}.cis_qtl.txt.gz"

if [[ ! -f "${CIS}" ]]; then
  log_message "ERROR: missing cis QTL results: ${CIS}"
  exit 1
fi

log_message "Summarizing meQTL QC for ${REGION}"
python "${H}/05_qc_summarize.py" \
  --region "${REGION}" \
  --cis-qtl "${CIS}"

conda deactivate
log_message "**** Job ends ****"
