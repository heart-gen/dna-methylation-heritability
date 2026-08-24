#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomicsguest
#SBATCH --job-name=cmb-meqtl
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-22
#SBATCH --time=06:00:00
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#
# One array task per autosome. meqtl_parameters.yml fixes chromosomes:
# autosomal, so there is no chrX task and none is silently skipped -- the array
# bound is the policy.

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

RUN_ID=${CMB_RUN_ID:?CMB_RUN_ID must be set}
CHROM=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID must be set}
RUN_DIR="${REPO_DIR}/05_cpg_meqtl_burden/_m/runs/${RUN_ID}"

log_job_info
require_file "${RUN_DIR}/results/tested-cpg-membership.tsv"

# tensorqtl lives in its own environment; it is not in the epigenomics env and
# must not be pulled in there, because it drags a pinned torch build.
V2_ENV_TENSORQTL="${V2_ENV_TENSORQTL:-$ENV_PATH/tensorqtl}"
require_file "$V2_ENV_TENSORQTL"

log_message "mapping cis-meQTL for chr${CHROM}"

conda run --no-capture-output -p "$V2_ENV_TENSORQTL" \
    python "${REPO_DIR}/05_cpg_meqtl_burden/_h/02_map_cpg_meqtl.py" \
    --run-id "$RUN_ID" \
    --chrom "$CHROM" \
    --threads "$V2_THREADS"

log_message "chr${CHROM} complete"
