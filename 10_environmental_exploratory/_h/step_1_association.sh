#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --job-name=env-assoc
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --array=1-22
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#
# One task per autosome. This stage reads only per-locus methylation phenotypes
# and a donor-level exposure table -- no genotype, no dosage matrix -- so it is
# far lighter than the SNP-consuming modules and the ceiling is set by the
# number of VMRs on the chromosome rather than by variant count.

# SLURM copies the batch script to /var/spool, so ${BASH_SOURCE[0]} does NOT
# resolve to _h/ at run time. V2_REPO_ROOT is exported by the submit driver;
# fall back to the submit directory for a hand-run sbatch.
V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

RUN_ID="${ENV_RUN_ID:?ENV_RUN_ID must be exported}"
CHROM="${SLURM_ARRAY_TASK_ID:?this script must run as a SLURM array over 1-22}"
# Analysis code runs from the run's immutable snapshot, not live _h/.
H_DIR="${V2_RUN_CODE:-${REPO_DIR}/10_environmental_exploratory/_h}"

log_job_info
log_message "chr${CHROM}, run ${RUN_ID}"

run_r "${H_DIR}/02_vmr_exposure_association.R" \
    --run-id "$RUN_ID" --chrom "$CHROM"
log_message "chr${CHROM} exposure association slice complete"
