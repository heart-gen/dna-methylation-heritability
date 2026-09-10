# 01b_estimation_cells/_m — generated output only

Nothing in this directory is written by hand. Every file here is produced by a
script in `01b_estimation_cells/_h/` from locked configuration and declared
inputs (AGENTS.md §5.2).

## What is tracked in Git

- This README.

Everything else — `runs/`, PLINK output, logs — is gitignored and lives on
Quest. A single cell carries roughly 56,000 output files (three PLINK files per
VMR plus per-task extraction logs), so the run tree is large.

Note the gitignore patterns for module output are `[0-9][0-9]*_*/_m/...`, not
`[0-9][0-9]_*/_m/...`. This module's directory name has a third character of
`b`, and the narrower pattern did not match it.

## Run directories

Output lands in `runs/{RUN_ID}/`, named
`estcell-{cell}-{region}-{YYYYMMDD}`, where `cell` is the estimation-cell token
`{catalog_cohort}.{estimation_group}` (for example
`estcell-all_individuals.EA-dlpfc-20260910`). Run directories are **immutable**:
`new_run()` refuses to reuse an existing one, and `close_run()` sets the
contents read-only. Never update a completed run in place; make a new one.

A cell run mirrors the Module 01 directory contract, because
`00_shared/locus_io.R` reads it: `vmr/` (donor list, `vmr.bed`,
`vmr_catalog.tsv`, per-VMR phenotypes), `covs/` (the catalog's covariates
restricted to this cell, plus `genotype_pcs.tsv`), and `plink_format/`
(per-VMR cis genotype for this cell's donors). It also carries `manifest.tsv`,
`task_reconciliation.tsv`, `extraction-log.tsv` and `output_checksums.tsv`.

The loci are copied verbatim from the upstream Module 01 run and checksum-
verified at close. A cell never re-derives VMRs, phenotypes or residuals — it
re-partitions donors over a fixed locus set.

## Regenerating a run

1. Check out the `git_commit` recorded in its `manifest.tsv`.
2. Confirm the `config_*_sha256` fields match the current `config/`.
3. Re-run `_h/submit_estimation_cell.sh` from `01b_estimation_cells/_m`.

A completed SLURM job is not proof of scientific validity. Check
`task_reconciliation.tsv`, then `_h/05_report_acceptance.R`, then the module
README's acceptance gate.
