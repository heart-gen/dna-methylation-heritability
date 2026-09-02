#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --job-name=lsp_qc_pred
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%A.out
#SBATCH --error=%x-%A.err
#
# Module 03 step 6: post-hoc QC on the predictability distribution.
#
# Diagnostic, NOT an acceptance gate -- it runs against runs that are already
# sealed and accepted, and it writes to _m/qc/{run_id}/, never into the
# immutable run directory (AGENTS.md 5.2).
#
# 32G / 4 cpus and 2h: the LD-pruning pass reads ~15.4M autosomal variants for
# the KING kinship matrix, which is the only expensive part. The annotation
# overlaps and the 2 GB mappability bigwig are cheap because only the VMR
# intervals are read.
#
# Env-driven with no positional arguments:
#   LSP_RUN_ID=lsp-AA-caudate-20260825 \
#   LSP_COMPARE_RUNS=lsp-AA-dlpfc-20260825,lsp-AA-hippocampus-20260825 \
#     sbatch _h/step_6_qc_predictability.sh
#
# LSP_COMPARE_RUNS is optional; without it the cross-region concordance check is
# skipped and only the annotation, relatedness and tail-ladder checks run.

# SLURM copies the batch script to /var/spool, so ${BASH_SOURCE[0]} does NOT
# resolve to _h/ at run time. V2_REPO_ROOT is exported by the submit driver;
# fall back to the submit directory for a hand-run sbatch.
V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

RUN_ID=${LSP_RUN_ID:?LSP_RUN_ID must be set}

# Unlike stages 1-5 this stage deliberately runs from LIVE _h/, not from the
# run's code snapshot. The QC is authored after the run sealed, so the snapshot
# does not contain this script -- and the question being asked is about the
# sealed outputs, not about reproducing the code that made them.
H_DIR="${REPO_DIR}/03_local_snp_prediction/_h"

# plink2 is not on the default PATH and is not in the conda envs.
export PLINK2="${PLINK2:-/projects/p32505/opt/bin/plink2}"

ARGS=(--run-id "$RUN_ID")
if [ -n "${LSP_COMPARE_RUNS:-}" ]; then
    ARGS+=(--compare-runs "$LSP_COMPARE_RUNS")
fi
if [ -n "${LSP_SKIP_RELATEDNESS:-}" ]; then
    ARGS+=(--skip-relatedness)
fi

log_job_info
run_r "${H_DIR}/06_qc_predictability_artifact.R" "${ARGS[@]}" ${LSP_EXTRA_ARGS:-}
