#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=24G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=ea_m3a_burden
#SBATCH --output=logs/ea_m3a_burden.%j.log

# EA-M3a VMR burden (does not overwrite EA M0 under _m/EA/).
# Submit from meqtl-validation/02_vmr_meqtl_burden/_m/:
#   sbatch --export=ALL,REGION=caudate ../_h/step_ea_m3a_burden.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h"
REGION="${REGION:-caudate}"
POPULATION=EA

CIS="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m/tensorqtl/${POPULATION}/M3a/qc/lead_snp_per_cpg.tsv.gz"
OUTDIR="${ROOT}/meqtl-validation/02_vmr_meqtl_burden/_m/EA_M3a/${REGION}"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

[[ -f "${CIS}" ]] || { log_message "ERROR: missing ${CIS}"; exit 1; }

log_message "Aggregate EA-M3a burden ${REGION}"
python3 "${H}/01_aggregate_vmr_burden.py" \
  --region "${REGION}" \
  --population "${POPULATION}" \
  --cis-qtl "${CIS}" \
  --outdir "${OUTDIR}" \
  --fdr 0.05

log_message "Fit EA-M3a burden models ${REGION}"
python3 "${H}/02_fit_burden_models.py" \
  --region "${REGION}" \
  --burden-tsv "${OUTDIR}/vmr_meqtl_burden.tsv.gz"

log_message "**** Job ends ****"
