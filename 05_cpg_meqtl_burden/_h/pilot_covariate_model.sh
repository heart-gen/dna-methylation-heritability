#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --job-name=cmb-cov-pilot
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Decision pilot for the 05 covariate model (executed 3 snpPCs vs locked M3a),
# one autosome, four arms. NOT a production run: no run ID, nothing under _m/,
# and deliberately not named step_*.sh so the submit driver cannot pick it up.
# Outputs go to $PILOT_OUT in scratch. Run `methpcs` once before this.
#
#   sbatch --chdir=$PILOT_OUT/slurm \
#       --export=ALL,REPO=<repo>,PILOT_OUT=<scratch dir>,CHROM=10 \
#       05_cpg_meqtl_burden/_h/pilot_covariate_model.sh
set -euo pipefail

REPO=${REPO:?set REPO to the worktree root}
OUT=${PILOT_OUT:?set PILOT_OUT to a scratch directory}
CHROM=${CHROM:-10}
ENV=/projects/p32505/opt/envs/genomics
PILOT="$REPO/05_cpg_meqtl_burden/_h/pilot_covariate_model.py"

# Same thread ceiling `00_shared/slurm.sh` exports for the production step.
# Without it torch opens one thread per core on the node while SLURM has
# allocated 8, and the permutation pass runs about 3x slower.
CPUS=${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS="$CPUS" OPENBLAS_NUM_THREADS="$CPUS" MKL_NUM_THREADS="$CPUS"

for ARM in ${PILOT_ARMS:-executed m0 m3a executed_keep m3a_keep}; do
    echo "=== $(date -Is) prepare $ARM chr$CHROM ==="
    conda run --no-capture-output -p "$ENV" python "$PILOT" prepare \
        --arm "$ARM" --chrom "$CHROM" --out "$OUT" --threads 8
    echo "=== $(date -Is) map $ARM chr$CHROM ==="
    conda run --no-capture-output -p "$ENV" python "$PILOT" map \
        --arm "$ARM" --chrom "$CHROM" --out "$OUT" --threads 8
done

echo "=== $(date -Is) summarize chr$CHROM ==="
conda run --no-capture-output -p "$ENV" python "$PILOT" summarize \
    --chrom "$CHROM" --out "$OUT"
echo "=== $(date -Is) done ==="
