#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=0-2
#SBATCH --job-name=predcmp_aggregate

set -euo pipefail

: "${PREDICTOR_COMPARISON_RUN_ROOT:?PREDICTOR_COMPARISON_RUN_ROOT must be set}"
: "${PREDICTOR_COMPARISON_CODE_ROOT:?PREDICTOR_COMPARISON_CODE_ROOT must be set}"
: "${PREDICTOR_COMPARISON_REPO_ROOT:?PREDICTOR_COMPARISON_REPO_ROOT must be set}"

REGIONS=(caudate dlpfc hippocampus)
REGION=${REGION:-${REGIONS[${SLURM_ARRAY_TASK_ID}]}}
ENV_PATH=${PREDICTOR_COMPARISON_ENV:-/projects/p32505/opt/envs/genomics}
SCRIPT_DIR=${PREDICTOR_COMPARISON_CODE_ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h
LEAD=${PREDICTOR_COMPARISON_REPO_ROOT}/meqtl-validation/01_cpg_meqtl_mapping/${REGION}/_m/tensorqtl/qc/lead_snp_per_cpg.tsv.gz
OUTDIR=${PREDICTOR_COMPARISON_RUN_ROOT}/regions/${REGION}

mkdir -p "${OUTDIR}"
"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/01_aggregate_vmr_burden.py" \
    --region "${REGION}" \
    --population AA \
    --cis-qtl "${LEAD}" \
    --fdr 0.05 \
    --outdir "${OUTDIR}"

