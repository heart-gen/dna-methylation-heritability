#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --job-name=cmb-subsample
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
# No --array here: the driver passes it. Task i maps autosome (i % 22) + 1 for
# draw block (i / 22), so one array covers every chromosome x draw block.
#SBATCH --time=04:00:00
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

# See step_2_map_meqtl.sh for why V2_REPO_ROOT is resolved this way and why the
# code runs from the run's snapshot.
V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

RUN_ID=${CMB_RUN_ID:?CMB_RUN_ID must be set}
B=${CMB_DRAWS_N:?CMB_DRAWS_N must be set}
PER_TASK=${CMB_DRAWS_PER_TASK:?CMB_DRAWS_PER_TASK must be set}
TASK=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID must be set}
H_DIR="${V2_RUN_CODE:-${REPO_DIR}/05_cpg_meqtl_burden/_h}"

CHROM=$(( TASK % 22 + 1 ))
BLOCK=$(( TASK / 22 ))
START=$(( BLOCK * PER_TASK + 1 ))
END=$(( START + PER_TASK - 1 ))
[ "$END" -gt "$B" ] && END=$B

log_job_info
log_message "chr${CHROM}, draws ${START}-${END}"
conda run --no-capture-output -p "${V2_ENV_PY}" \
    python "${H_DIR}/07_subsample_map.py" \
    --run-id "$RUN_ID" --chrom "$CHROM" \
    --draw-start "$START" --draw-end "$END"
log_message "chr${CHROM} draws ${START}-${END} complete"
