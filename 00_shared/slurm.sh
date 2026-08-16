#!/bin/bash
# Shared SLURM helpers for v2 revision launchers.
#
# Source this at the top of every step_*.sh:
#
#     source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"
#
# Defect V12: legacy launchers used `module load plink` and another user's
# absolute paths. v2 never touches the module system -- the genomics conda env
# ships plink 1.9, which cannot read .pgen/.pvar, so plink2 comes from the shared
# opt tree.

set -euo pipefail

# ---------------------------------------------------------------- repo root
# Resolved the same way here::here() and 00_shared/load.R do, so shell and R
# agree on what "the repository" means.
REPO_DIR="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$REPO_DIR" != "/" ] && [ ! -d "$REPO_DIR/.git" ]; do
    REPO_DIR=$(dirname "$REPO_DIR")
done
if [ ! -d "$REPO_DIR/.git" ]; then
    echo "ERROR: could not locate repository root from ${SLURM_SUBMIT_DIR:-$PWD}" >&2
    exit 1
fi
export V2_REPO_ROOT="$REPO_DIR"

# ------------------------------------------------------------- environments
ENV_PATH="${ENV_PATH:-/projects/p32505/opt/envs}"
V2_ENV_R="${V2_ENV_R:-$ENV_PATH/epigenomics}"
V2_ENV_CALIBRATION="${V2_ENV_CALIBRATION:-$ENV_PATH/calibrated-local-h2}"
PLINK2="${PLINK2:-/projects/p32505/opt/bin/plink2}"

# Never inherit a half-configured module environment.
if command -v module >/dev/null 2>&1; then
    module purge >/dev/null 2>&1 || true
fi

# ------------------------------------------------------------------ threads
# Every threaded library under R sizes its pool from the NODE's core count, not
# from the cgroup SLURM gave us. On a Quest compute node that means data.table
# opened 128 threads inside a 1-CPU allocation (detectCores 256, getDTthreads
# 128) for the whole first full-scale DLPFC run. It did not cause a failure, but
# it oversubscribes the allocation and makes runtimes unreproducible.
#
# The fix has to be applied on BOTH sides. These variables cover the libraries
# that read the environment (OpenMP, and the BLAS implementations); data.table
# ignores them once its pool is built, so 00_shared/threads.R calls
# setDTthreads() as well. Neither alone is sufficient.
#
# V2_CPUS is the allocation. V2_THREADS is the ceiling anything may open: 2x the
# allocation, which allows modest hyperthread oversubscription without the
# runaway. Pass V2_THREADS to any tool with its own thread flag (plink2 --threads).
V2_CPUS="${SLURM_CPUS_PER_TASK:-${V2_CPUS:-1}}"
V2_THREADS=$(( V2_CPUS * 2 ))
export V2_CPUS V2_THREADS
export OMP_NUM_THREADS="$V2_CPUS"
export OPENBLAS_NUM_THREADS="$V2_CPUS"
export MKL_NUM_THREADS="$V2_CPUS"
export R_DATATABLE_NUM_THREADS="$V2_CPUS"

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_job_info() {
    echo "**** QUEST info ****"
    echo "User:       ${USER:-unknown}"
    echo "Job id:     ${SLURM_JOB_ID:-none}"
    echo "Job name:   ${SLURM_JOB_NAME:-none}"
    echo "Array task: ${SLURM_ARRAY_TASK_ID:-none}"
    echo "Hostname:   ${HOSTNAME:-unknown}"
    echo "Repo:       ${REPO_DIR}"
    echo "CPUs:       ${V2_CPUS} (thread ceiling ${V2_THREADS})"
    echo "Commit:     $(git -C "$REPO_DIR" rev-parse --short HEAD 2>/dev/null || echo unknown)"
    echo "********************"
}

# Verify a required executable or file exists before burning array tasks on it.
require_file() {
    if [ ! -e "$1" ]; then
        echo "ERROR: required path missing: $1" >&2
        exit 1
    fi
}

require_exec() {
    if [ ! -x "$1" ]; then
        echo "ERROR: required executable not found: $1" >&2
        exit 1
    fi
}

# Run an R script in the v2 R environment, from the repository root.
run_r() {
    local script="$1"; shift
    require_file "$script"
    conda run --no-capture-output -p "$V2_ENV_R" Rscript "$script" "$@"
}

# Per-chromosome length lookup.
#
# Defect V10: legacy step_5.sh ran
#     CHR_SIZE=$(grep "^chr1[[:space:]]" $CHR_FILE | cut -f2)
# once and used chr1's length as the bound for EVERY chromosome. chr1 is the
# longest, so the check never fired on smaller chromosomes.
chrom_size() {
    local chr="$1"
    local chr_file="${2:-/projects/b1213/resources/genomes/human/gencode-v47/fasta/chromosome_sizes.txt}"
    require_file "$chr_file"
    local size
    size=$(awk -v c="$chr" '$1 == c || $1 == "chr"c {print $2; exit}' "$chr_file")
    if [ -z "$size" ]; then
        echo "ERROR: no size entry for chromosome '$chr' in $chr_file" >&2
        exit 1
    fi
    echo "$size"
}
