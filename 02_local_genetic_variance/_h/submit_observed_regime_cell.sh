#!/bin/bash

## Submit the observed-regime diagnostic grid for one cohort-by-region cell
## that has no completed production run.
##
## Usage: submit_observed_regime_cell.sh COHORT REGION VMR_RUN_ID [DATE]
##
## Chain: locus-geometry array -> combine -> scenario manifest -> scenario
## array (fixed 1844 chunks) -> summary. The geometry scan exists because the
## grid stratifies loci on num_snps and p_eff, which are pure genotype
## geometry: computing them needs no phenotype, elastic net or BSLMM, so a
## cell can be characterised without first running a full production pass.
##
## This grid does not change the frozen model, its calibration, or its gates.

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: $0 COHORT REGION VMR_RUN_ID [DATE]" >&2
    exit 2
fi
COHORT=$1
REGION=$2
VMR_RUN_ID=$3
DATE=${4:-$(date +%Y%m%d)}

H_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MODULE_DIR=$(cd "${H_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
R_BIN=${ENV_PATH}/bin/Rscript
ACCOUNT=${SBATCH_ACCOUNT:-p32505}
PARTITION=${LGV_PARTITION:-short}
MAX_CONCURRENT=${LGV_MAX_CONCURRENT:-100}
LOCI_PER_CHUNK=${LGV_LOCI_PER_CHUNK:-25}

GEOMETRY_RUN_ID="lgv-geometry-${COHORT}-${REGION}-${DATE}"
REGIME_RUN_ID="lgv-observed-regime-${COHORT}-${REGION}-${DATE}"
RUNS_ROOT="${MODULE_DIR}/_m/runs"
REGIME_DIR="${RUNS_ROOT}/${REGIME_RUN_ID}"

GEOMETRY_DIR=$("${R_BIN}" "${H_DIR}/11_prepare_locus_geometry.R" \
    --run-id="${GEOMETRY_RUN_ID}" --cohort="${COHORT}" --region="${REGION}" \
    --vmr-run-id="${VMR_RUN_ID}" --vmrs-per-chunk="${LOCI_PER_CHUNK}")
GEOMETRY_DIR=$(echo "${GEOMETRY_DIR}" | tail -1 | tr -d '[:space:]')
N_GEOMETRY=$(awk -F '\t' 'NR > 1 {if ($1 > n) n=$1} END {print n+0}' \
    "${GEOMETRY_DIR}/config/chunk-manifest.tsv")

## The scenario count is fixed by the locked config (loci_per_stratum x strata
## x architectures x PVE levels x replicates), so the scenario array can be
## sized now rather than after the manifest job runs.
N_SCENARIO=$("${R_BIN}" -e '
a <- commandArgs(TRUE)
s <- read.delim(a[[1]], stringsAsFactors = FALSE)
v <- function(k) s$value[s$setting == k]
loci <- as.integer(v("loci_per_stratum")) * as.integer(v("num_snps_strata")) *
    as.integer(v("p_eff_strata"))
n <- loci * length(strsplit(v("architectures"), ",")[[1]]) *
    length(strsplit(v("h2_values"), ",")[[1]]) *
    as.integer(v("replicates_per_cell"))
cat(ceiling(n / as.integer(v("scenarios_per_chunk"))))
' "${MODULE_DIR}/config/observed-regime-20260822.tsv")

EXPORTS="ALL,LGV_GEOMETRY_DIR=${GEOMETRY_DIR},LGV_RUN_DIR=${REGIME_DIR}"
EXPORTS="${EXPORTS},LGV_H_DIR=${H_DIR},CAL_H2_ENV=${ENV_PATH}"
EXPORTS="${EXPORTS},LGV_REGIME_RUN_ID=${REGIME_RUN_ID},LGV_COHORT=${COHORT}"
EXPORTS="${EXPORTS},LGV_REGION=${REGION},LGV_VMR_RUN_ID=${VMR_RUN_ID}"

sub() {
    local dependency=$1 step=$2
    shift 2
    sbatch --parsable --account="${ACCOUNT}" --partition="${PARTITION}" \
        ${dependency:+--dependency="${dependency}"} \
        --job-name="lgv_${COHORT}_${REGION}_${step%.sh}" \
        --output="${GEOMETRY_DIR}/logs/%x.%j.log" \
        --export="${EXPORTS}" "$@" "${H_DIR}/${step}" | cut -d';' -f1
}

GEO_JOB=$(sbatch --parsable --account="${ACCOUNT}" --partition="${PARTITION}" \
    --array="1-${N_GEOMETRY}%${MAX_CONCURRENT}" \
    --job-name="lgv_${COHORT}_${REGION}_geometry" \
    --output="${GEOMETRY_DIR}/logs/%x.%A_%a.log" \
    --export="${EXPORTS}" "${H_DIR}/step_09_locus_geometry.sh" | cut -d';' -f1)

## afterany: a cancelled or failed geometry array must still be reconciled.
COMB_JOB=$(sub "afterany:${GEO_JOB}" step_10_combine_geometry.sh)
MAN_JOB=$(sub "afterok:${COMB_JOB}" step_11_regime_manifest.sh)
SCEN_JOB=$(sbatch --parsable --account="${ACCOUNT}" --partition="${PARTITION}" \
    --dependency="afterok:${MAN_JOB}" \
    --array="1-${N_SCENARIO}%${MAX_CONCURRENT}" \
    --job-name="lgv_${COHORT}_${REGION}_regime" \
    --output="${GEOMETRY_DIR}/logs/%x.%A_%a.log" \
    --export="${EXPORTS}" \
    "${H_DIR}/step_07_observed_regime_scenarios.sh" | cut -d';' -f1)
SUM_JOB=$(sub "afterany:${SCEN_JOB}" step_08_summarize_observed_regime.sh)

printf 'stage\tstep_script\tjob_id\n' > "${GEOMETRY_DIR}/submitted-jobs.tsv"
printf '%s\t%s\t%s\n' \
    geometry  step_09_locus_geometry.sh            "${GEO_JOB}" \
    combine   step_10_combine_geometry.sh          "${COMB_JOB}" \
    manifest  step_11_regime_manifest.sh           "${MAN_JOB}" \
    scenarios step_07_observed_regime_scenarios.sh "${SCEN_JOB}" \
    summary   step_08_summarize_observed_regime.sh "${SUM_JOB}" \
    >> "${GEOMETRY_DIR}/submitted-jobs.tsv"

cat <<MSG
Submitted observed-regime grid for ${COHORT}/${REGION}
  geometry run : ${GEOMETRY_DIR} (${N_GEOMETRY} chunks)
  regime run   : ${REGIME_DIR} (${N_SCENARIO} chunks)
  jobs         : geometry=${GEO_JOB} combine=${COMB_JOB} manifest=${MAN_JOB}
                 scenarios=${SCEN_JOB} summary=${SUM_JOB}
MSG
