#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=cal_h2_vmr
#SBATCH --output=calibrated-simulation-analysis/_m/logs/%x.%A_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=04:00:00

set -euo pipefail

if [[ -n "${CAL_H2_SCRIPT_DIR:-}" ]]; then
    SCRIPT_DIR=$(readlink -f "${CAL_H2_SCRIPT_DIR}")
elif [[ -n "${SLURM_SUBMIT_DIR:-}" &&
        -d "${SLURM_SUBMIT_DIR}/calibrated-simulation-analysis/_h" ]]; then
    SCRIPT_DIR=$(readlink -f \
        "${SLURM_SUBMIT_DIR}/calibrated-simulation-analysis/_h")
else
    SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
fi
ANALYSIS_DIR=${CAL_H2_ANALYSIS_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}
if [[ ! -f "${SCRIPT_DIR}/04_estimate_observed_vmr.R" ]]; then
    echo "Analysis script is missing from SCRIPT_DIR: ${SCRIPT_DIR}" >&2
    exit 1
fi
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
REPO_ROOT=${CAL_H2_REPO_ROOT:-$(cd "${ANALYSIS_DIR}/.." && pwd)}
CALIBRATION_MODEL=${CAL_H2_CALIBRATION_MODEL:?CAL_H2_CALIBRATION_MODEL must be set}
OUTPUT_ROOT=${CAL_H2_OBSERVED_OUTPUT_ROOT:?CAL_H2_OBSERVED_OUTPUT_ROOT must be set}
PLINK_ROOT=${CAL_H2_PLINK_ROOT:-/projects/b1213/users/alexis/projects/dna-methylation-heritability/vmr-analysis/all_individuals}
RECOVERED_PLINK_ROOT=${CAL_H2_RECOVERED_PLINK_ROOT:-}
PGEN_ROOT=${CAL_H2_PGEN_ROOT:-/projects/b1213/users/alexis/projects/dna-methylation-heritability/inputs/genotypes/all_individuals}
PHENOTYPE_ROOT=${CAL_H2_PHENOTYPE_ROOT:-/projects/b1213/users/alexis/projects/dna-methylation-heritability/vmr-analysis/all_individuals}
WRITE_DIAGNOSTICS=${CAL_H2_WRITE_DIAGNOSTICS:-FALSE}
: "${REGION:?REGION must be set}"
: "${POPULATION:?POPULATION must be set}"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

ARRAY_INDEX=${SLURM_ARRAY_TASK_ID}
TASK_ID=${ARRAY_INDEX}
if [[ -n "${CAL_H2_TASK_MANIFEST:-}" ]]; then
    if [[ ! -s "${CAL_H2_TASK_MANIFEST}" ]]; then
        echo "Task manifest is missing: ${CAL_H2_TASK_MANIFEST}" >&2
        exit 1
    fi
    TASK_ID=$(awk -F '\t' -v row="$((ARRAY_INDEX + 1))" \
        'NR == row {print $1}' "${CAL_H2_TASK_MANIFEST}")
    if [[ ! "${TASK_ID}" =~ ^[1-9][0-9]*$ ]]; then
        echo "Could not resolve task ID for recovery array index ${ARRAY_INDEX}" >&2
        exit 1
    fi
fi

# A missing per-VMR PLINK file is an upstream computational gap, not a locus-QC
# outcome. During an explicitly configured recovery run, reconstruct that exact
# prespecified cis window from the immutable cohort PGEN. A genuine zero-variant
# extraction is marked for the R adapter to record as a QC failure.
if [[ -n "${RECOVERED_PLINK_ROOT}" ]]; then
    VMR_FILE=${REPO_ROOT}/vmr-analysis/all_individuals/${REGION,,}/_m/vmr.bed
    VMR_RECORD=$(awk -v row="${TASK_ID}" 'NR == row {print $1, $2, $3}' "${VMR_FILE}")
    read -r CHROMOSOME START END <<< "${VMR_RECORD}"
    if [[ -z "${CHROMOSOME:-}" || ! "${START:-}" =~ ^[0-9]+$ ||
          ! "${END:-}" =~ ^[0-9]+$ ]]; then
        echo "Could not resolve VMR coordinates for task ${TASK_ID}" >&2
        exit 1
    fi
    CHROMOSOME_LABEL=${CHROMOSOME#chr}
    STEM=${START}_${END}
    UPSTREAM_PREFIX=${PLINK_ROOT}/${REGION,,}/_m/plink_format/chr_${CHROMOSOME_LABEL}/TOPMed_LIBD-${POPULATION}.${STEM}
    RECOVERY_DIR=${RECOVERED_PLINK_ROOT}/${REGION,,}/_m/plink_format/chr_${CHROMOSOME_LABEL}
    RECOVERY_PREFIX=${RECOVERY_DIR}/TOPMed_LIBD-${POPULATION}.${STEM}
    UPSTREAM_NO_SNPS=FALSE
    if [[ -s "${UPSTREAM_PREFIX}.log" ]] &&
        grep -Fq 'No variants remaining after main filters' "${UPSTREAM_PREFIX}.log"; then
        UPSTREAM_NO_SNPS=TRUE
    fi
    if [[ ! -s "${UPSTREAM_PREFIX}.bed" && ! -s "${RECOVERY_PREFIX}.bed" &&
          ! -e "${RECOVERY_PREFIX}.no-snps" &&
          "${UPSTREAM_NO_SNPS}" != TRUE ]]; then
        mkdir -p "${RECOVERY_DIR}"
        WINDOW_START=$((START - 500000))
        if (( WINDOW_START < 1 )); then WINDOW_START=1; fi
        WINDOW_END=$((END + 500000))
        KEEP_FILE=${PHENOTYPE_ROOT}/${REGION,,}/_m/samples-${POPULATION}.txt
        if [[ ! -s "${PGEN_ROOT}/TOPMed_LIBD.pgen" ||
              ! -s "${PGEN_ROOT}/TOPMed_LIBD.pvar" ||
              ! -s "${PGEN_ROOT}/TOPMed_LIBD.psam" || ! -s "${KEEP_FILE}" ]]; then
            echo "Cannot reconstruct missing PLINK window: cohort PGEN or keep file is missing" >&2
            exit 1
        fi
        # V12: plink2 from the shared opt tree, never the module system. The
        # genomics conda env ships plink 1.9, which cannot read .pgen/.pvar.
        PLINK2="${PLINK2:-/projects/p32505/opt/bin/plink2}"
        if [ ! -x "$PLINK2" ]; then
            echo "ERROR: plink2 not found at $PLINK2" >&2
            exit 1
        fi
        set +e
        "$PLINK2" --pfile "${PGEN_ROOT}/TOPMed_LIBD" \
            --threads "${SLURM_CPUS_PER_TASK:-1}" \
            --chr "${CHROMOSOME_LABEL}" \
            --from-bp "${WINDOW_START}" \
            --to-bp "${WINDOW_END}" \
            --make-bed \
            --keep "${KEEP_FILE}" \
            --no-parents --no-sex --no-pheno \
            --out "${RECOVERY_PREFIX}"
        plink_status=$?
        set -e
        if (( plink_status != 0 )); then
            if [[ -s "${RECOVERY_PREFIX}.log" ]] &&
                grep -Fq 'No variants remaining after main filters' "${RECOVERY_PREFIX}.log"; then
                : > "${RECOVERY_PREFIX}.no-snps"
            else
                echo "PLINK recovery extraction failed for task ${TASK_ID}" >&2
                exit "${plink_status}"
            fi
        fi
    fi
fi

set +e
"${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/04_estimate_observed_vmr.R" \
    --region="${REGION}" \
    --population="${POPULATION}" \
    --task-id="${TASK_ID}" \
    --repo-root="${REPO_ROOT}" \
    --plink-root="${PLINK_ROOT}" \
    --recovered-plink-root="${RECOVERED_PLINK_ROOT}" \
    --phenotype-root="${PHENOTYPE_ROOT}" \
    --calibration-model="${CALIBRATION_MODEL}" \
    --output-root="${OUTPUT_ROOT}" \
    --write-diagnostics="${WRITE_DIAGNOSTICS}"
status=$?
set -e
if (( status != 0 )); then
    failure_dir=${OUTPUT_ROOT}/${REGION,,}/${POPULATION}/failures
    mkdir -p "${failure_dir}"
    printf 'task_id\tregion\tpopulation\texit_status\tlog_file\n%s\t%s\t%s\t%s\t%s\n' \
        "${TASK_ID}" "${REGION,,}" "${POPULATION}" "${status}" \
        "${SLURM_JOB_NAME:-cal_h2_vmr}.${SLURM_ARRAY_JOB_ID:-NA}_${ARRAY_INDEX}.log" \
        > "${failure_dir}/vmr-$(printf '%07d' "${TASK_ID}").tsv"
    echo "Recorded failed VMR task ${TASK_ID} with exit status ${status}" >&2
fi
exit 0
