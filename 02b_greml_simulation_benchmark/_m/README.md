# 02b_greml_simulation_benchmark/_m — generated output only

Nothing in this directory is written by hand. Every file here is produced by a
script in `02b_greml_simulation_benchmark/_h/` from locked configuration and declared inputs.

## What is tracked, and where the rest lives

Files are tiered by `supplementary_data/_h/build_release.py`, and every file of
every accepted run is listed with its tier in `supplementary_data/release_manifest.tsv`:

- **git**: this README; `combined/` (the cross-region deliverables built from
  accepted runs, except `-UNACCEPTED` outputs and bulk external reference data);
  small provenance files a script asserts against at runtime; and, from the
  accepted runs, every file the manuscript is written from (figure source
  tables, Supplementary Data tables, decision and gate tables, `manifest.tsv`)
  that is 15 MB or smaller. Those run files are added with `git add -f`, since
  `runs/` itself stays gitignored.
- **Git LFS**: the same class of manuscript file above 15 MB.
- **Zenodo** (one zip per module, DOI 10.5281/zenodo.20547606): reproducibility
  extras no manuscript text uses, plus a mirror of the two tiers above.
- **Not distributed**: individual-level data (methylation matrices, genotypes,
  covariates; dbGaP phs000979) and bulk intermediates (per-task shards, logs,
  full nominal output), which are regenerable from `_h/` and stay on Quest.

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

Each run also holds `inputs/` (arm 1 input checksums), `config/` (task
manifests, scenario grid, REML modes), `status/` (one row per unit, which is
what reconciliation counts), `results/` (per-task fits), `summary/` (estimates,
recovery metrics, Spearman, gate criteria) and `figures/`. Arm 1 keeps its
combined LD-score table and SNP group lists in `work/`; the GRMs built from
them and all per-task scratch are removed before the seal and the removal is
recorded in `manifest.tsv` as `removed_before_seal`.

## Combined output

`combined/` holds `10_collate.R`'s tables over the accepted runs:
`greml-benchmark-{recovery-metrics,spearman,decisions}.tsv`. A table built with
`--allow-unaccepted TRUE` carries `-UNACCEPTED` in its name and
`citable = FALSE`.

## Regenerating a run

1. Check out the `git_commit` recorded in its `manifest.tsv`.
2. Confirm the `config_*_sha256` fields match the current `config/`.
3. Re-run `02b_greml_simulation_benchmark/_h/submit_greml_benchmark.sh ar1`, or
   `... observed AA <region>`.

A completed SLURM job is not proof of scientific validity. Check
`summary/gate-criteria.tsv` and the module README's acceptance gate.
