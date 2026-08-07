#!/bin/bash
# One-time move of LIBD eQTL mapping outputs/scripts from Phase 7 into this module.
# Safe to run only AFTER the in-flight TensorQTL job writing under 08 has completed.
set -euo pipefail

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
SRC="${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_m/level3/libd_eqtl/caudate"
DST="${ROOT}/meqtl-validation/09_libd_eqtl_mapping/_m/caudate"
H08="${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_h"
H09="${ROOT}/meqtl-validation/09_libd_eqtl_mapping/_h"
CIS_SRC="${SRC}/genes/tensorqtl/libd_aa_caudate_genes_standard.cis_qtl.txt.gz"

echo "$(date '+%Y-%m-%d %H:%M:%S') - Reorganizing LIBD eQTL into 09_libd_eqtl_mapping"

if [[ ! -f "${CIS_SRC}" ]]; then
  echo "ERROR: TensorQTL cis output not found: ${CIS_SRC}"
  echo "Wait for TensorQTL to finish before reorganizing."
  exit 1
fi

if [[ -e "${DST}" && ! -L "${DST}" ]]; then
  echo "ERROR: destination already exists: ${DST}"
  exit 1
fi

mkdir -p "$(dirname "${DST}")"
mv "${SRC}" "${DST}"
echo "Moved ${SRC} -> ${DST}"

# Compatibility symlink so any old Level-3 paths still resolve during transition
mkdir -p "${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_m/level3/libd_eqtl"
ln -sfn "${DST}" "${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_m/level3/libd_eqtl/caudate"
echo "Symlink: 08/.../level3/libd_eqtl/caudate -> ${DST}"

# Remove obsolete mapping scripts from 08 (canonical copies live in 09)
for f in 14_prepare_libd_aa_expression.R 15_make_libd_aa_covariates.R 16_make_tensorqtl_bed.py \
         step_6b_libd_eqtl_prep.sh step_6c_libd_eqtl_tensorqtl.sh; do
  if [[ -f "${H08}/${f}" ]]; then
    rm -f "${H08}/${f}"
    echo "Removed obsolete 08 script: ${f}"
  fi
done

# Record provenance
printf 'moved_from\tmoved_to\tcis_qtl\tdate\n%s\t%s\t%s\t%s\n' \
  "${SRC}" "${DST}" "${DST}/genes/tensorqtl/libd_aa_caudate_genes_standard.cis_qtl.txt.gz" \
  "$(date -Iseconds)" \
  > "${ROOT}/meqtl-validation/09_libd_eqtl_mapping/_m/REORGANIZED_FROM_08.tsv"

echo "$(date '+%Y-%m-%d %H:%M:%S') - Done"
