#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --job-name=sldsc_ld
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#
# Stage 05: annotate one chromosome and compute its LD scores.
#
# One array task per autosome. Each task builds its own .annot.gz and then the
# matching .l2.ldscore.gz, so a failed chromosome is retried alone rather than
# rerunning all 22.

source "$(dirname "${BASH_SOURCE[0]}")/../../00_shared/slurm.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_ID="${SLDSC_RUN_ID:?SLDSC_RUN_ID must be exported}"
CHR="${SLURM_ARRAY_TASK_ID:?this step runs as a SLURM array}"
RUN_DIR="${REPO_DIR}/06_partitioned_heritability/_m/runs/${RUN_ID}"

log_job_info
log_message "chromosome ${CHR}, run ${RUN_ID}"

PY_ENV="/projects/p32505/opt/envs/genomics"

# pybedtools writes large intermediates; keep them off /var/tmp on compute nodes.
export TMPDIR="${TMPDIR:-/projects/b1213/tmp}"
mkdir -p "$TMPDIR"

LDSC_DIR=$(conda run --no-capture-output -p "$PY_ENV" python -c \
    "import yaml,sys; print(yaml.safe_load(open('${REPO_DIR}/config/partitioned_heritability.yml'))['ldsc_dir'])")
require_file "$LDSC_DIR/ldsc.py"

BIM_PREFIX=$(conda run --no-capture-output -p "$PY_ENV" python -c \
    "import yaml; c=yaml.safe_load(open('${REPO_DIR}/config/partitioned_heritability.yml')); r=c['ld_references'][c['ld_reference_arm']]; print(r['bim_dir']+'/'+r['bim_prefix'])")
PRINT_SNPS=$(conda run --no-capture-output -p "$PY_ENV" python -c \
    "import yaml; c=yaml.safe_load(open('${REPO_DIR}/config/partitioned_heritability.yml')); r=c['ld_references'][c['ld_reference_arm']]; print(r['hapmap3_print_snps'])")
LD_WIND=$(conda run --no-capture-output -p "$PY_ENV" python -c \
    "import yaml; print(yaml.safe_load(open('${REPO_DIR}/config/partitioned_heritability.yml'))['ld_wind_cm'])")

require_file "${BIM_PREFIX}${CHR}.bim"
require_file "$PRINT_SNPS"

log_message "building annotation for chr${CHR}"
conda run --no-capture-output -p "$PY_ENV" python \
    "${SCRIPT_DIR}/03_make_annot.py" --run-id "$RUN_ID" --chrom "$CHR"

ANNOT="${RUN_DIR}/ldscores/annot.${CHR}.annot.gz"
require_file "$ANNOT"

log_message "computing LD scores for chr${CHR}"
conda run --no-capture-output -p "$PY_ENV" python \
    "${SCRIPT_DIR}/ldsc_wrapper.py" "$LDSC_DIR" ldsc.py \
    --l2 \
    --bfile "${BIM_PREFIX}${CHR}" \
    --ld-wind-cm "$LD_WIND" \
    --annot "$ANNOT" \
    --thin-annot \
    --print-snps "$PRINT_SNPS" \
    --out "${RUN_DIR}/ldscores/annot.${CHR}"

require_file "${RUN_DIR}/ldscores/annot.${CHR}.l2.ldscore.gz"
log_message "chr${CHR} complete"
