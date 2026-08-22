# 02_local_genetic_variance/_m — generated output only

Nothing in this directory is written by hand. Every file here is produced by a
script in `02_local_genetic_variance/_h/` from locked configuration and declared inputs
(AGENTS.md §5.2).

## What is tracked in Git

- This README.
- Small provenance files a script asserts against at runtime, where they exist
  (for example `02_local_genetic_variance/_m/calibration_frozen/`).

Everything else — `runs/`, matrices, PLINK output, figures, logs — is
gitignored and lives on Quest.

## Run directories

Output lands in `runs/{RUN_ID}/`, named
`{module}-{cohort}-{region}-{YYYYMMDD}`. Run directories are **immutable**:
`new_run()` refuses to reuse an existing one, and `close_run()` sets the
contents read-only. Never update a completed run in place; make a new one.

Each run carries `manifest.tsv` (git commit, config checksums, input checksums,
upstream run IDs, `vmr_set_id`, ordered donor checksum, sample counts, seeds,
software environment, SLURM job IDs),
`results/combined/task-reconciliation.tsv`, and `output_checksums.tsv`.

Observed runs additionally carry `config/task-manifest.tsv` and
`config/chunk-manifest.tsv`, which map every VMR to exactly one initial array
chunk. Stage 01 writes one terminal row per VMR under `results/task_rows/`.
Stage 02 reconciles those rows against the frozen task universe and converts a
missing scheduler output into an explicit computational failure. Production
acceptance therefore requires zero computational failures and exact task
reconciliation.

An active observed run must also carry one combined
`local-genetic-control-*-vmrs.tsv` produced by
`_h/04_derive_local_snp_contribution_score.R`. This is the downstream contract:
eligible VMRs have a within-cell rank score and standardized score, while every
excluded VMR retains an explicit reason. Raw joint estimates may remain in
audit tables but are not accepted as exact locus-level PVE.

## Regenerating a run

1. Check out the `git_commit` recorded in its `manifest.tsv`.
2. Confirm the `config_*_sha256` fields match the current `config/`.
3. Re-run the active workflow with a new immutable run ID:
   `_h/submit_observed_local_control.sh`.

A completed SLURM job is not proof of scientific validity. Check
`results/combined/task-reconciliation.tsv`,
`results/combined/observed-score-decision.tsv`, and the module README's
acceptance gate.

The estimator-development, boundary-remapping, validation, recovery, and old
observed-EN submitters are archived under
`_h/archive/2026-08-21_estimator_development/`. They are forensic records, not
supported execution paths.
