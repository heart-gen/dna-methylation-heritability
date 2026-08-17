#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=cal_h2_audit_sum
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

set -euo pipefail
SCRIPT_DIR=${CAL_H2_SCRIPT_DIR:?CAL_H2_SCRIPT_DIR must be set}
RUN_ROOT=${CAL_H2_AUDIT_RUN_ROOT:?CAL_H2_AUDIT_RUN_ROOT must be set}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/09_summarize_high_signal_audit.R" \
    --manifest="${RUN_ROOT}/config/high-signal-audit-manifest.tsv" \
    --summaries="${RUN_ROOT}/summaries" \
    --output-dir="${RUN_ROOT}/summary"
