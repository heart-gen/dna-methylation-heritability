#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --array=0-65
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=meqtl_cpg_cov
#SBATCH --output=logs/cpg_coverage.%A_%a.log

# 3 regions x 22 autosomes = 66 tasks (0-65)
# Submit from meqtl-validation/01_cpg_meqtl_mapping/_m:
#   mkdir -p logs && sbatch ../_h/step_6_coverage.sh

set -euo pipefail
ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_h"
cd "${ROOT}/meqtl-validation/01_cpg_meqtl_mapping/_m"
mkdir -p logs

REGIONS=(caudate dlpfc hippocampus)
CHROMS=({1..22})
TASK=${SLURM_ARRAY_TASK_ID}
REGION_IDX=$((TASK / 22))
CHROM_IDX=$((TASK % 22))
REGION="${REGIONS[$REGION_IDX]}"
CHROM="${CHROMS[$CHROM_IDX]}"
POPULATION="${POPULATION:-AA}"

echo "$(date '+%Y-%m-%d %H:%M:%S') - coverage extract ${REGION} chr${CHROM} (${POPULATION})"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/epigenomics

Rscript "${H}/06_extract_cpg_coverage.R" \
  --region "${REGION}" \
  --chrom "${CHROM}" \
  --population "${POPULATION}" \
  --min-coverage 5

echo "$(date '+%Y-%m-%d %H:%M:%S') - done ${REGION} chr${CHROM} (${POPULATION})"
