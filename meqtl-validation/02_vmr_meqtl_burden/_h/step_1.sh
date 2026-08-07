#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --array=0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_vmr_aggregate
#SBATCH --output=logs/vmr_aggregate.%A_%a.log

# Submit from meqtl-validation/02_vmr_meqtl_burden/_m/:
#   mkdir -p logs && sbatch ../_h/step_1.sh

set -euo pipefail

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
mkdir -p logs

REGIONS=(caudate dlpfc hippocampus)
if [[ -z "${REGION:-}" ]]; then
  REGION="${REGIONS[${SLURM_ARRAY_TASK_ID}]}"
fi
POPULATION="${POPULATION:-AA}"

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h"
if [[ "${POPULATION}" == "AA" ]]; then
  CIS="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m/tensorqtl/qc/lead_snp_per_cpg.tsv.gz"
else
  CIS="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m/tensorqtl/${POPULATION}/qc/lead_snp_per_cpg.tsv.gz"
fi

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Aggregating CpG meQTL burden for ${REGION} (${POPULATION})"
python3 "${H}/01_aggregate_vmr_burden.py" \
  --region "${REGION}" \
  --population "${POPULATION}" \
  --cis-qtl "${CIS}" \
  --fdr 0.05

log_message "**** Job ends ****"
