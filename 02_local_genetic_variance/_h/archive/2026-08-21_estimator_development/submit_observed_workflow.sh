#!/bin/bash

set -euo pipefail

echo "This legacy single-region entry point is disabled." >&2
echo "Use submit_observed.sh RUN_ID <AA|all_individuals> <region>." >&2
echo "The replacement enforces calibration acceptance, snapshots provenance, and submits sequential VMR chunks." >&2
exit 2
