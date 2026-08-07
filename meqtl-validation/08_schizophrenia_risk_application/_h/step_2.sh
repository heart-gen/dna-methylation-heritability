#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=scz_risk_meqtl
#SBATCH --output=logs/scz_risk_meqtl.%j.log

# Targeted risk-variant–CpG meQTL (Analysis 3) + architecture (4) + tx (6 light).
# Submit from _m/ after step_1:
#   sbatch --dependency=afterok:<step1> --export=ALL,REGION=caudate ../_h/step_2.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_h"
REGION="${REGION:-caudate}"
POPULATION="${POPULATION:-AA}"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Analysis 3: targeted risk-variant CpG meQTL (${REGION}, ${POPULATION})"
python3 "${H}/03_test_risk_variant_cpg_meqtl.py" \
  --region "${REGION}" \
  --population "${POPULATION}"

log_message "Analysis 4: predictability architecture"
python3 "${H}/04_architecture_predictability.py" --region "${REGION}"

log_message "Analysis 6 (light): transcriptional coupling"
python3 "${H}/05_tx_integration.py" --region "${REGION}"

log_message "**** Job ends ****"
