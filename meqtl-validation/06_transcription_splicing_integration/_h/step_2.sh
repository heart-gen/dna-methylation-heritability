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
#SBATCH --job-name=meqtl_tx_psi
#SBATCH --output=logs/tx_psi.%A_%a.log

# Submit from meqtl-validation/06_transcription_splicing_integration/_m/:
#   sbatch ../_h/step_2.sh

set -euo pipefail

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

REGIONS=(caudate dlpfc hippocampus)
if [[ -z "${REGION:-}" ]]; then
  REGION="${REGIONS[${SLURM_ARRAY_TASK_ID}]}"
fi

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/06_transcription_splicing_integration/_h"
BURDEN="${ROOT}/meqtl-validation/02_vmr_meqtl_burden/_m/${REGION}/vmr_meqtl_burden.tsv.gz"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "meQTL × PSI enrichment for ${REGION}"
python3 "${H}/01_meqtl_tx_enrichment.py" \
  --region "${REGION}" \
  --modality psi \
  --burden-tsv "${BURDEN}"

# optional ABC (small N; exploratory)
if [[ -f "${ROOT}/heritability/elastic_net_model/all_individuals/tissue_comparison/regulatory_context/_m/${REGION}/AA/expression/abc/architecture_model_input.tsv" ]]; then
  log_message "meQTL × ABC expression enrichment for ${REGION}"
  python3 "${H}/01_meqtl_tx_enrichment.py" \
    --region "${REGION}" \
    --modality expression_abc \
    --burden-tsv "${BURDEN}" || log_message "WARNING: ABC enrichment failed (often tiny N)"
fi

log_message "**** Job ends ****"
