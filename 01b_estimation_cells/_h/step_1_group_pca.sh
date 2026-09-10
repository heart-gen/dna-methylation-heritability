#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=estcell_pca
#SBATCH --output=logs/estcell_pca.%j.log
#
# 01b step 1: genotype PCs computed WITHIN this estimation cell's donors.
#
# Why not reuse inputs/genotypes/all_individuals/TOPMed_LIBD.eigenvec: those
# eigenvectors were computed on the pooled sample, where the leading PCs
# separate the donor groups from each other. Restricted to one group they are
# close to constant and carry little within-group information, so they cannot
# control the structure that remains INSIDE the group. PI decision 2026-09-10
# (config/covariates.yml: estimation_cells.genotype_pcs).
#
# The PCA is run on the same pooled pgen the cis windows are extracted from, so
# the variants are identical; only the donors are restricted.
#
# Usage:
#   cd 01b_estimation_cells/_m && mkdir -p logs
#   RUN_ID=<id> sbatch ../_h/step_1_group_pca.sh

_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$_ROOT" != "/" ] && [ ! -d "$_ROOT/.git" ]; do _ROOT=$(dirname "$_ROOT"); done
source "$_ROOT/00_shared/slurm.sh"

: "${RUN_ID:?set RUN_ID}"

RUN_DIR="$REPO_DIR/01b_estimation_cells/_m/runs/$RUN_ID"
require_file "$RUN_DIR/manifest.tsv"
require_file "$RUN_DIR/vmr/donors_plink.txt"
require_exec "$PLINK2"

manifest_value() {
    awk -F'\t' -v f="$1" '$1 == f { print $2; found = 1 } END { if (!found) exit 1 }' \
        "$RUN_DIR/manifest.tsv"
}

PFILE=$(manifest_value pgen_prefix)
N_PC=$(manifest_value genotype_pcs | awk -F',' '{print NF}')
: "${PFILE:?manifest carries no pgen_prefix}"
require_file "${PFILE}.pgen"

N_DONORS=$(wc -l < "$RUN_DIR/vmr/donors_plink.txt")

log_job_info
log_message "within-group PCA: $N_DONORS donors, $N_PC PCs, pfile $PFILE"

# plink2 needs more variance-standardized markers than PCs; with n donors it can
# return at most n-1 components. Fail here rather than silently returning fewer
# columns than the covariate model expects.
if [ "$N_DONORS" -le "$N_PC" ]; then
    echo "ERROR: $N_DONORS donors cannot support $N_PC PCs" >&2
    exit 1
fi

WORK="$RUN_DIR/covs/pca"
mkdir -p "$WORK"

# LD pruning first. Without it the leading PCs track a handful of long-range LD
# blocks (the MHC, the chr8 and chr17 inversions) rather than genome-wide
# structure, which is the classic way a PCA covariate adjusts for the wrong
# thing. --maf is applied WITHIN this cell's donors, which is the point.
"$PLINK2" --pfile "$PFILE" \
          --keep "$RUN_DIR/vmr/donors_plink.txt" \
          --autosome \
          --snps-only just-acgt \
          --maf 0.05 \
          --geno 0.05 \
          --hwe 1e-6 \
          --indep-pairwise 200kb 1 0.1 \
          --threads "$V2_THREADS" \
          --out "$WORK/prune"

require_file "$WORK/prune.prune.in"
log_message "pruned to $(wc -l < "$WORK/prune.prune.in") variants"

"$PLINK2" --pfile "$PFILE" \
          --keep "$RUN_DIR/vmr/donors_plink.txt" \
          --extract "$WORK/prune.prune.in" \
          --pca "$N_PC" approx \
          --threads "$V2_THREADS" \
          --out "$WORK/group"

require_file "$WORK/group.eigenvec"

# Normalize to the schema load_observed_locus() reads: FID IID snpPC1..snpPCk.
# plink2 writes '#FID IID PC1 PC2 ...' (or '#IID ...' for a single-column-ID
# fileset), so the header is rewritten explicitly rather than assumed.
run_r "$REPO_DIR/01b_estimation_cells/_h/01_write_genotype_pcs.R" \
    --run-id "$RUN_ID"

require_file "$RUN_DIR/covs/genotype_pcs.tsv"
log_message "wrote $RUN_DIR/covs/genotype_pcs.tsv"
log_message "**** Job ends ****"
