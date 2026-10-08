#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --qos=buyin
#SBATCH --job-name=sldsc_ext_ld
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err
#
# External-annotation benchmark, stage 11: for ONE chromosome, build every
# staged annotation's thin-annot column and its LD scores, with the same
# reference bim, ld_wind_cm and print-snps as 05_compute_ldscores.sh.
#
# Annotations are computed one file each (not one multi-column file) because
# each is regressed ALONE on top of baselineLD in stage 12, and LDSC reads every
# column of a --ref-ld-chr file.

V2_REPO_ROOT="${V2_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [ "$V2_REPO_ROOT" != "/" ] && [ ! -d "$V2_REPO_ROOT/.git" ]; do
    V2_REPO_ROOT=$(dirname "$V2_REPO_ROOT")
done
source "${V2_REPO_ROOT}/00_shared/slurm.sh"

RUN_ID="${SLDSC_RUN_ID:?SLDSC_RUN_ID must be exported}"
SCRIPT_DIR="${V2_RUN_CODE:-${REPO_DIR}/06_partitioned_heritability/_h}"
CHR="${SLURM_ARRAY_TASK_ID:?this step runs as a SLURM array}"
RUN_DIR="${REPO_DIR}/06_partitioned_heritability/_m/runs/${RUN_ID}"
CFG="${RUN_DIR}/code/config/partitioned_heritability.yml"
require_file "$CFG"
LIST="${RUN_DIR}/annotation/annotation-list.tsv"
require_file "$LIST"

log_job_info
PY_ENV="/projects/p32505/opt/envs/genomics"
PY="conda run --no-capture-output -p ${PY_ENV} python"
export TMPDIR="${TMPDIR:-/projects/b1213/tmp}"
mkdir -p "$TMPDIR"

read -r LDSC_DIR BIM_PREFIX PRINT_SNPS LD_WIND < <($PY -c "
import yaml; c=yaml.safe_load(open('${CFG}')); r=c['ld_references'][c['ld_reference_arm']]
print(c['ldsc_dir'], r['bim_dir']+'/'+r['bim_prefix'], r['hapmap3_print_snps'], c['ld_wind_cm'])")
require_file "$LDSC_DIR/ldsc.py"
require_file "${BIM_PREFIX}${CHR}.bim"

ANNOTS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next} $h["status"]=="ready"{print $h["annotation"]}' "$LIST")
for A in $ANNOTS; do
    OUT="${RUN_DIR}/ldscores/${A}/annot.${CHR}"
    log_message "${A} chr${CHR}: annotation"
    $PY "${SCRIPT_DIR}/11a_external_make_annot.py" --run-id "$RUN_ID" --chrom "$CHR" --annotation "$A"
    require_file "${OUT}.annot.gz"
    log_message "${A} chr${CHR}: LD scores"
    $PY "${SCRIPT_DIR}/ldsc_wrapper.py" "$LDSC_DIR" ldsc.py \
        --l2 --bfile "${BIM_PREFIX}${CHR}" --ld-wind-cm "$LD_WIND" \
        --annot "${OUT}.annot.gz" --thin-annot --print-snps "$PRINT_SNPS" \
        --out "$OUT"
    require_file "${OUT}.l2.ldscore.gz"
    $PY "${SCRIPT_DIR}/11a_external_make_annot.py" --run-id "$RUN_ID" --chrom "$CHR" \
        --annotation "$A" --relabel
done
log_message "chr${CHR} complete"
