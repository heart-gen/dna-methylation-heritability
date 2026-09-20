# 09b_aging_application/_m — generated output only

Nothing in this directory is written by hand. Every file here is produced by a
script in `09b_aging_application/_h/` from locked configuration and declared inputs
(AGENTS.md §5.2).

## What is tracked in Git

- This README.
- Small provenance files a script asserts against at runtime, where they exist
  (for example `02_local_genetic_variance/_m/calibration_frozen/`).
- `combined/`: the cross-region summary tables, built from accepted runs. These
  are the module's deliverable -- the decision rows, the concordance tables, and
  the per-VMR result tables a journal would ask for as Supplementary Data -- so
  the numbers behind a claim are readable from a clone without re-running SLURM.
  Two things under `combined/` are not tracked: outputs of a stage whose config
  is not PI-locked, which carry `-UNACCEPTED` in the filename and may not be
  cited, and bulk external reference data such as the public GWAS summary
  statistics under `09_schizophrenia_risk_application`.

Everything else — `runs/`, matrices, PLINK output, figures, logs — is
gitignored and lives on Quest.

## Run directories

Output lands in `runs/{RUN_ID}/`, named
`{module}-{cohort}-{region}-{YYYYMMDD}`. Run directories are **immutable**:
`new_run()` refuses to reuse an existing one, and `close_run()` sets the
contents read-only. Never update a completed run in place; make a new one.

Each run carries `manifest.tsv` (git commit, config checksums, upstream run
IDs, `vmr_set_id`, ordered donor checksum, sample counts, bootstrap count,
software environment, SLURM job IDs) and `output_checksums.tsv`. This module has
no task array, so there is no `task_reconciliation.tsv`; non-finite VMRs are
listed in `results/excluded-vmrs.tsv` instead.

`runs/{RUN_ID}/checkpoint/` holds the donor x VMR methylation matrix, the
design matrices, the axis designs and the bootstrap draws. They stay with the
run because stage 05 bootstraps against them and the draws are half of the
headline variance.

## Cross-region output

`combined/` holds stage 05's tables, `aging-*-{cohort}.tsv`. A table built with
`--allow-unlocked` carries `-UNACCEPTED` in its name and `citable = FALSE`.

## Regenerating a run

1. Check out the `git_commit` recorded in its `manifest.tsv`.
2. Confirm the `config_*_sha256` fields match the current `config/`.
3. Re-run `09b_aging_application/_h/submit_aging.sh AA <region>`.

A completed SLURM job is not proof of scientific validity. Check
`results/gate-checks.tsv` and the module README's acceptance gate.
