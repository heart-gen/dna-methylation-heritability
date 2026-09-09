#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --job-name=scz-coloc
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --array=1-22
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#
# One task per autosome: prepare the harmonised regions, then run coloc.abf over
# them. The two stages are chained inside one task rather than across two array
# jobs because the regions are large intermediate files local to the chromosome,
# and a dependency between two 22-way arrays would serialise on the slowest
# chromosome twice.
#
# Preparation runs in the genomics env (pyarrow); coloc runs in the coloc env.
#
# Sizing is measured, not guessed. chr22 (3 loci, 13 QTL resources) peaked at
# ~5 GB and took ~11 min for preparation alone. Loci per chromosome vary about
# tenfold (3 on chr22, 35 on chr1), and preparation re-reads each resource once
# per locus chunk, so the headroom here covers chr1. SCZ_COLOC_LOCUS_CHUNK
# tunes that trade-off: smaller chunks use less memory and more I/O.
export SCZ_COLOC_LOCUS_CHUNK="${SCZ_COLOC_LOCUS_CHUNK:-5}"

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

run_py "${H_DIR}/08_prepare_coloc_regions.py" --run-id "$RUN_ID" --chrom "$CHROM"
log_message "chr${CHROM} regions prepared"

run_r_coloc "${H_DIR}/09_run_coloc.R" --run-id "$RUN_ID" --chrom "$CHROM"
log_message "chr${CHROM} coloc.abf complete"
