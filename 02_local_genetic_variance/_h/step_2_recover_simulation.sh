#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=cal_h2_sim_recovery
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=04:00:00

set -euo pipefail
SCRIPT_DIR=${CAL_H2_SCRIPT_DIR:?CAL_H2_SCRIPT_DIR must be set}
MANIFEST=${SCENARIO_MANIFEST:?SCENARIO_MANIFEST must be set}
RECOVERY_MANIFEST=${CAL_H2_RECOVERY_MANIFEST:?CAL_H2_RECOVERY_MANIFEST must be set}
OUTPUT_ROOT=${SIMULATION_OUTPUT_ROOT:?SIMULATION_OUTPUT_ROOT must be set}
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
CHUNK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID must be set}
mapfile -t SCENARIO_TASK_IDS < <(awk -F '\t' -v chunk="${CHUNK_ID}" \
    'NR > 1 && $1 == chunk {print $2}' "${RECOVERY_MANIFEST}")
if (( ${#SCENARIO_TASK_IDS[@]} == 0 )); then
    echo "Chunk ${CHUNK_ID} has no scenario IDs" >&2
    exit 1
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
for SCENARIO_TASK_ID in "${SCENARIO_TASK_IDS[@]}"; do
    if [[ ! "${SCENARIO_TASK_ID}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Chunk ${CHUNK_ID} contains an invalid scenario ID" >&2
        exit 1
    fi
    "${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/01_simulate_and_crossfit.R" \
        --manifest="${MANIFEST}" \
        --task-id="${SCENARIO_TASK_ID}" \
        --output-root="${OUTPUT_ROOT}"
done
