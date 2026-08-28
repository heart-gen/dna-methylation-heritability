#!/bin/bash
## 10 / Table 1 and cohort QC (PI decision D3).
## Usage: bash step_2_table1_qc.sh RUN_ID
set -euo pipefail
source "$(dirname "$0")/../../00_shared/slurm.sh"

RUN_ID="${1:?usage: step_2_table1_qc.sh RUN_ID}"

log_message "[10] Table 1 and ancestry panel"
run_r "10_integrated_manuscript_outputs/_h/04_table1_cohort.R" --run-id "$RUN_ID"

log_message "[10] cross-region sample-integrity screen"
run_r "10_integrated_manuscript_outputs/_h/05_qc_sample_integrity.R" --run-id "$RUN_ID"

log_message "[10] done: $RUN_ID"
