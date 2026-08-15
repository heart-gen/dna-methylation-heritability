#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --array=0-2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_cov_aggregate
#SBATCH --output=logs/cpg_coverage_aggregate.%A_%a.log

# Aggregate chromosome-level coverage after the full step_6 array succeeds.

set -euo pipefail
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h"
cd "${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_m"
mkdir -p logs

REGIONS=(caudate dlpfc hippocampus)
REGION="${REGION:-${REGIONS[${SLURM_ARRAY_TASK_ID}]}}"
POPULATION="${POPULATION:-AA}"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics
python3 "${H}/07_aggregate_vmr_coverage.py" \
  --region "${REGION}" \
  --population "${POPULATION}"
