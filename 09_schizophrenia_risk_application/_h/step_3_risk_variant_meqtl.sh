#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --job-name=scz-rvmeqtl
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --array=1-22
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#
# One task per autosome. The memory ceiling is set by Module 05's nominal
# parquet for the chromosome (~0.5-0.7 GB on disk, several times that in
# memory); the stage pushes both filters into the parquet scan so the whole
# file is never materialised.

# SLURM copies the batch script to /var/spool, so ${BASH_SOURCE[0]} does NOT
# resolve to _h/ at run time. V2_REPO_ROOT is exported by the submit driver;
# fall back to the submit directory for a hand-run sbatch.
V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

RUN_ID="${SCZ_RUN_ID:?SCZ_RUN_ID must be exported}"
CHROM="${SLURM_ARRAY_TASK_ID:?this script must run as a SLURM array over 1-22}"
# Analysis code runs from the run's immutable snapshot, not live _h/.
H_DIR="${V2_RUN_CODE:-${REPO_DIR}/09_schizophrenia_risk_application/_h}"

log_job_info
log_message "chr${CHROM}, run ${RUN_ID}"

run_py "${H_DIR}/03_risk_variant_cpg_meqtl.py" \
    --run-id "$RUN_ID" --chrom "$CHROM" --threads "${V2_CPUS}"
log_message "chr${CHROM} risk-variant meQTL slice complete"
