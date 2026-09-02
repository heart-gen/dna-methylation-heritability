#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --job-name=tsc-assoc
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# One job per modality. The PSI assay is ~691k events x 487 libraries; running
# it on the submit host OOM-killed the stage, hence the large --mem here and the
# feature subsetting inside 02_run_local_associations.R.

# SLURM copies the batch script to /var/spool, so ${BASH_SOURCE[0]} does NOT
# resolve to _h/ at run time. V2_REPO_ROOT is exported by the submit driver;
# fall back to the submit directory for a hand-run sbatch.
V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

RUN_ID="${TSC_RUN_ID:?TSC_RUN_ID must be exported}"
MODALITY="${TSC_MODALITY:?TSC_MODALITY must be exported}"
# Analysis code runs from the run's immutable snapshot, not live _h/.
H_DIR="${V2_RUN_CODE:-${REPO_DIR}/07_transcription_splicing_coupling/_h}"

log_job_info
log_message "modality ${MODALITY}, run ${RUN_ID}"

run_r "${H_DIR}/02_run_local_associations.R" --run-id "$RUN_ID" --modality "$MODALITY"
log_message "modality ${MODALITY} complete"
