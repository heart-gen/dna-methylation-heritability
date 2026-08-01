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
#SBATCH --job-name=meqtl_calib
#SBATCH --output=logs/cpg_calibration.%A_%a.log

# Submit from meqtl-validation/01_cpg_meqtl_mapping/_m:
#   mkdir -p logs && sbatch ../_h/step_7_calibration.sh
# Significance: Storey qval <= FDR (default 0.05)

set -euo pipefail
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h"
cd "${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_m"
mkdir -p logs

REGIONS=(caudate dlpfc hippocampus)
REGION="${REGION:-${REGIONS[${SLURM_ARRAY_TASK_ID}]}}"
FDR="${FDR:-0.05}"
CIS="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m/tensorqtl/cpg_meqtl_${REGION}.cis_qtl.txt.gz"

echo "$(date '+%Y-%m-%d %H:%M:%S') - calibration ${REGION}"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

python3 "${H}/08_calibration_plots.py" \
  --region "${REGION}" \
  --cis-qtl "${CIS}" \
  --fdr "${FDR}"

echo "$(date '+%Y-%m-%d %H:%M:%S') - done ${REGION}"
