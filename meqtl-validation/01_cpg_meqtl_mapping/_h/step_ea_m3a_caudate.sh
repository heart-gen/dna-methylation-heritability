#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=ea_m3a_latent
#SBATCH --output=logs/ea_m3a_latent.%j.log

# Build EA-specific methPC1–5 and covariates_M3a for caudate.
# Submit from meqtl-validation/01_cpg_meqtl_mapping/_m/:
#   sbatch ../_h/step_ea_m3a_caudate.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h"
REGION="${REGION:-caudate}"
POPULATION=EA

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

log_message "Estimate EA methPCs for ${REGION}"
python3 "${H}/09_estimate_latent_factors.py" --region "${REGION}" --population "${POPULATION}"

log_message "Build EA covariate models (incl. M3a)"
python3 "${H}/10_prepare_covariate_models.py" --region "${REGION}" --population "${POPULATION}"

# Promote EA M3a covariates to a dedicated path for TensorQTL
PREP="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m/prepared/${POPULATION}"
/bin/cp -f "${PREP}/covariates_M3a.txt" "${PREP}/covariates_M3a_primary.txt"
log_message "Wrote ${PREP}/covariates_M3a_primary.txt"
head -1 "${PREP}/covariates_M3a_primary.txt"

log_message "**** Job ends ****"
