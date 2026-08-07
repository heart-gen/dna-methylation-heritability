#!/bin/bash

set -euo pipefail

echo "This legacy single-region entry point is disabled." >&2
echo "Use submit_observed_calibrated_workflow.sh CALIBRATION_RUN_ID OBSERVED_RUN_ID [AA]." >&2
echo "The replacement enforces calibration acceptance, snapshots provenance, submits all three regions, and runs completeness/QC aggregation." >&2
exit 2
