#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=96G
#SBATCH --mail-type=FAIL
#SBATCH --array=1-66
#SBATCH --job-name=vmr_universe
#SBATCH --output=logs/universe.%A_%a.log
#SBATCH --error=logs/universe.%A_%a.log

## Decompose the caudate VMR/CpG universe advantage: sample size or coverage
## depth? See ../_h/01_universe_decomposition.R for the full rationale.
##
## Array spans 22 autosomes x 3 regions. Task IDs 1-22 caudate, 23-44 dlpfc,
## 45-66 hippocampus. X/Y are excluded -- VMR calling uses autosomes and the
## CT-SNP masks are autosome-only.

set -euo pipefail

ENV_PATH="/projects/p32505/opt/envs/epigenomics"

TASK=$((SLURM_ARRAY_TASK_ID - 1))
REGION_IDX=$((TASK / 22))
CHR=$(((TASK % 22) + 1))

case "$REGION_IDX" in
    0) REGION="caudate" ;;
    1) REGION="dlpfc" ;;
    2) REGION="hippocampus" ;;
    *) echo "Invalid task id: $SLURM_ARRAY_TASK_ID" >&2; exit 1 ;;
esac

echo "**** Job starts $(date '+%Y-%m-%d %H:%M:%S') ****"
echo "Job id: ${SLURM_JOB_ID}  Task: ${SLURM_ARRAY_TASK_ID}"
echo "Region: ${REGION}  Chromosome: ${CHR}"

conda run --no-capture-output -p "$ENV_PATH" \
    Rscript ../_h/01_universe_decomposition.R "$CHR" "$REGION"

echo "**** Job ends $(date '+%Y-%m-%d %H:%M:%S') ****"
