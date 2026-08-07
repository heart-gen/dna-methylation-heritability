#!/bin/bash

# Submit the reproducible DNAm deconvolution DAG from any working directory.
set -euo pipefail

ROOT="/projects/b1213/users/kynon/projects/dna-methylation-heritability"
H="${ROOT}/inputs/cell_proportions/_h"
M="${ROOT}/inputs/cell_proportions/_m"
ENV_PREFIX="${CELL_DECONV_ENV:-${M}/conda_env}"
if [[ ! -x "${ENV_PREFIX}/bin/Rscript" ]]; then
  echo "Missing ${ENV_PREFIX}. From ${M}, run:" >&2
  echo "  conda env create --prefix ./conda_env -f ../../../conda_environments/cell_deconvolution.yml" >&2
  exit 1
fi
mkdir -p "${M}/logs"
cd "${M}"

job_id() { printf '%s' "${1%%;*}"; }
REF_JOB="$(job_id "$(sbatch --parsable "${H}/step_3_reference.sh")")"
EXTRACT_JOB="$(job_id "$(sbatch --parsable --dependency="afterok:${REF_JOB}" "${H}/step_4_extract.sh")")"
DECONV_JOB="$(job_id "$(sbatch --parsable --dependency="afterok:${EXTRACT_JOB}" "${H}/step_5_deconvolve.sh")")"
VALIDATE_JOB="$(job_id "$(sbatch --parsable --dependency="afterok:${DECONV_JOB}" "${H}/step_6_validate.sh")")"

{
  printf 'stage\tjob_id\tdependency\n'
  printf 'reference\t%s\tnone\n' "${REF_JOB}"
  printf 'extract\t%s\tafterok:%s\n' "${EXTRACT_JOB}" "${REF_JOB}"
  printf 'deconvolve\t%s\tafterok:%s\n' "${DECONV_JOB}" "${EXTRACT_JOB}"
  printf 'validate\t%s\tafterok:%s\n' "${VALIDATE_JOB}" "${DECONV_JOB}"
} > "${M}/dnam_scmd_submission_manifest.tsv"

printf 'Submitted DNAm scMD workflow: reference=%s extract=%s deconvolve=%s validate=%s\n' \
  "${REF_JOB}" "${EXTRACT_JOB}" "${DECONV_JOB}" "${VALIDATE_JOB}"
