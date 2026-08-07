#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --job-name=scz_l3_eqtl
#SBATCH --output=logs/scz_l3_eqtl.%j.log

# Level 3: targeted risk-variant eQTL + summary (consumes 09_libd_eqtl_mapping).
# Genome-wide LIBD eQTL prep/map is NOT run here — use:
#   meqtl-validation/09_libd_eqtl_mapping/_h/step_1.sh + step_2.sh
#
# Submit from 08 _m/ after step_6a and after 09 TensorQTL:
#   sbatch ../_h/step_6b_level3_eqtl_tests.sh

set -euo pipefail
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
mkdir -p logs

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_h"
CIS09="${ROOT}/meqtl-validation/09_libd_eqtl_mapping/_m/caudate/genes/tensorqtl/libd_aa_caudate_genes_standard.cis_qtl.txt.gz"
CIS08="${ROOT}/meqtl-validation/08_schizophrenia_risk_application/_m/level3/libd_eqtl/caudate/genes/tensorqtl/libd_aa_caudate_genes_standard.cis_qtl.txt.gz"

source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

CIS=""
if [[ -f "${CIS09}" ]]; then
  CIS="${CIS09}"
elif [[ -f "${CIS08}" ]]; then
  CIS="${CIS08}"
fi

CIS_ARG=()
if [[ -n "${CIS}" ]]; then
  log_message "Using cis-QTL: ${CIS}"
  CIS_ARG=(--cis-qtl "${CIS}")
else
  log_message "WARNING: cis_qtl not found yet; running targeted tests without eGene annotation"
fi

python3 "${H}/17_libd_risk_variant_eqtl.py" "${CIS_ARG[@]}"
python3 "${H}/18_summarize_level3.py"

log_message "**** Job ends ****"
