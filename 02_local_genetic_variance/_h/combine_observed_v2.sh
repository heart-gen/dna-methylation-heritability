#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=cal_h2_combine_v2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#
# 02_local_genetic_variance: combine and QC the v2 observed runs.
#
# The migrated step_6_combine_observed.sh handles exactly one observed run and
# is given its paths by the legacy submitter. v2 produces one run per arm x
# region cell, so this driver runs 07_combine_observed.R once per cell -- which
# is where the per-cell QC gate lives -- and then merges the per-cell tables
# into one module-level result keyed by vmr_set_id.
#
# Run directly or via sbatch, from the module's _m/ directory:
#   ../_h/combine_observed_v2.sh lgv-AA-caudate-20260816 [more run ids...]
#
# It combines whatever results exist, so it is useful on a partially finished
# run: a cell with pending tasks gets `complete FALSE` in its QC row and the
# whole driver exits non-zero. Read a non-zero exit as "not production yet",
# and the QC table for which cell and why.

set -euo pipefail

(( $# >= 1 )) || { echo "Usage: $0 OBSERVED_RUN_ID [OBSERVED_RUN_ID...]" >&2; exit 1; }

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS_DIR=$(cd "${SCRIPT_DIR}/.." && pwd)
REPO_ROOT=$(cd "${ANALYSIS_DIR}/.." && pwd)
ENV_PATH=${CAL_H2_ENV:-/projects/p32505/opt/envs/calibrated-local-h2}
COMBINED=${ANALYSIS_DIR}/_m/combined
mkdir -p "${COMBINED}"

# config/thresholds.yml gates.max_outside_calibration_domain, converted to the
# within-domain floor 07_combine_observed.R actually checks. Read from the
# config so the gate cannot drift from the number the PI locked.
MAX_OUTSIDE=$(awk '/^gates:/{g=1} g && /max_outside_calibration_domain:/{print $2; exit}' \
    "${REPO_ROOT}/config/thresholds.yml")
: "${MAX_OUTSIDE:?could not read gates.max_outside_calibration_domain}"
MIN_WITHIN=$(awk -v m="${MAX_OUTSIDE}" 'BEGIN{printf "%.4f", 1 - m}')

# fail_on_qc=FALSE here so that ONE failing cell does not hide the QC tables for
# the other five. Every cell's verdict is collected and the exit status is set
# from all of them at the end, which is the same fail-closed behaviour with a
# more useful report.
status=0
for run_id in "$@"; do
    ROOT=${ANALYSIS_DIR}/_m/observed-runs/${run_id}
    [[ -d "$ROOT" ]] || { echo "No such observed run: $ROOT" >&2; exit 1; }

    echo "=== ${run_id}"
    "${ENV_PATH}/bin/Rscript" "${SCRIPT_DIR}/07_combine_observed.R" \
        --input="${ROOT}/results" \
        --expected="${ROOT}/config/expected-tasks.tsv" \
        --output-dir="${ROOT}/results/combined" \
        --min-within-domain-rate="${MIN_WITHIN}" \
        --fail-on-qc=FALSE
done

# One module-level table across every cell. The per-cell files already carry
# region, population, upstream_vmr_run_id and vmr_set_id on every row, so the
# merge is a concatenation and nothing has to be inferred from a filename.
set +e
"${ENV_PATH}/bin/Rscript" - "$COMBINED" "$ANALYSIS_DIR" "$@" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
combined_dir <- args[[1L]]
analysis_dir <- args[[2L]]
run_ids <- args[-c(1L, 2L)]

read_one <- function(path) {
    if (!file.exists(path)) return(NULL)
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}
gather <- function(pattern) {
    tables <- lapply(run_ids, function(run_id) {
        dir <- file.path(analysis_dir, "_m", "observed-runs", run_id,
                         "results", "combined")
        hits <- list.files(dir, pattern = pattern, full.names = TRUE)
        if (!length(hits)) return(NULL)
        do.call(rbind, lapply(hits, read_one))
    })
    tables <- Filter(Negate(is.null), tables)
    if (!length(tables)) return(NULL)
    columns <- Reduce(union, lapply(tables, names))
    tables <- lapply(tables, function(d) {
        for (missing in setdiff(columns, names(d))) d[[missing]] <- NA
        d[, columns, drop = FALSE]
    })
    do.call(rbind, tables)
}

estimates <- gather("^calibrated-local-h2-.*-vmrs\\.tsv$")
qc <- gather("^observed-run-qc\\.tsv$")
if (is.null(estimates)) stop("No per-cell estimate tables were found")
if (is.null(qc)) stop("No per-cell QC tables were found")

## A locus is identified by its VMR catalog and its interval, never by task_id
## alone: task_id is a row number in one run's vmr.bed and means nothing across
## cells. Guard against a merge that silently collides two arms.
key <- paste(estimates$vmr_set_id, estimates$vmr_id)
if (anyDuplicated(key)) {
    stop("Duplicate (vmr_set_id, vmr_id) across cells: the merge is not 1:1")
}

dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
write.table(estimates, file.path(combined_dir, "calibrated-local-h2-all-cells.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(qc, file.path(combined_dir, "observed-run-qc-all-cells.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat("\n=== per-cell QC\n")
show <- c("region", "population", "vmr_set_id", "expected_tasks", "analyzed_tasks",
          "qc_failed_tasks", "computational_failed_tasks", "within_domain_rate",
          "complete", "overall_qc_pass")
print(qc[, intersect(show, names(qc))], row.names = FALSE)

cat("\n=== h2_en_calibrated by cell\n")
by_cell <- do.call(rbind, lapply(split(estimates, list(estimates$population, estimates$region), drop = TRUE), function(d) {
    data.frame(cell = paste(d$population[[1L]], d$region[[1L]]),
               n = nrow(d),
               median = round(median(d$h2_en_calibrated), 4),
               q25 = round(quantile(d$h2_en_calibrated, 0.25), 4),
               q75 = round(quantile(d$h2_en_calibrated, 0.75), 4),
               boundary_hits = sum(d$h2_upper_boundary_hit),
               stringsAsFactors = FALSE)
}))
print(by_cell, row.names = FALSE)

if (!all(qc$overall_qc_pass)) {
    cat("\nQC FAILED for:",
        paste(qc$region[!qc$overall_qc_pass], qc$population[!qc$overall_qc_pass],
              collapse = ", "), "\n")
    quit(save = "no", status = 2L)
}
RSCRIPT
status=$?
set -e

if (( status != 0 )); then
    echo "Combine completed but the QC gate did not pass." >&2
fi
exit "$status"
