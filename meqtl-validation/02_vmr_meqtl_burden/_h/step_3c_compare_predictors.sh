#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --array=0-2
#SBATCH --job-name=predcmp_models

set -euo pipefail

: "${PREDICTOR_COMPARISON_RUN_ROOT:?PREDICTOR_COMPARISON_RUN_ROOT must be set}"
: "${PREDICTOR_COMPARISON_CODE_ROOT:?PREDICTOR_COMPARISON_CODE_ROOT must be set}"
: "${PREDICTOR_COMPARISON_REPO_ROOT:?PREDICTOR_COMPARISON_REPO_ROOT must be set}"

REGIONS=(caudate dlpfc hippocampus)
REGION=${REGION:-${REGIONS[${SLURM_ARRAY_TASK_ID}]}}
ENV_PATH=${PREDICTOR_COMPARISON_ENV:-/projects/p32505/opt/envs/genomics}
SCRIPT_DIR=${PREDICTOR_COMPARISON_CODE_ROOT}/meqtl-validation/02_vmr_meqtl_burden/_h
OUTDIR=${PREDICTOR_COMPARISON_RUN_ROOT}/regions/${REGION}
BURDEN=${OUTDIR}/vmr_meqtl_burden.tsv.gz
LEGACY=${PREDICTOR_COMPARISON_REPO_ROOT}/heritability/elastic_net_model/all_individuals/${REGION}/_m/${REGION}_summary_elastic-net_AA.tsv
CONFIG=${PREDICTOR_COMPARISON_RUN_ROOT}/config/predictor_comparison.yml
CALIBRATED=${PREDICTOR_COMPARISON_RUN_ROOT}/config/calibrated-local-h2-AA-vmrs.tsv
CALIBRATED_QC=${PREDICTOR_COMPARISON_RUN_ROOT}/config/observed-run-qc.tsv

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/02_fit_burden_models.py" \
    --region "${REGION}" \
    --burden-tsv "${BURDEN}" \
    --seed 20260808 \
    --require-complete-tech-join

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/04_build_predictor_bridge.py" \
    --region "${REGION}" \
    --burden "${BURDEN}" \
    --legacy "${LEGACY}" \
    --calibrated "${CALIBRATED}" \
    --calibrated-qc "${CALIBRATED_QC}" \
    --config "${CONFIG}" \
    --output-dir "${OUTDIR}"

"${ENV_PATH}/bin/python" "${SCRIPT_DIR}/05_fit_predictor_comparison.py" \
    --region "${REGION}" \
    --bridge "${OUTDIR}/predictor_bridge.tsv.gz" \
    --config "${CONFIG}" \
    --output-dir "${OUTDIR}"

