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

## Known provenance defect in the nine existing runs (recorded 2026-09-23)

**Step 1 of "Regenerating a run" below does not work for six of the nine runs.**
This is a provenance defect, not a data defect, and because `_m/` is immutable it
cannot be repaired in place — it is recorded here so that a reader who tries to
regenerate a run meets it before wasting an afternoon on it.

All nine runs record `git_dirty = true`. Six — the three
`estcell-all_individuals.AA-*-20260910` and three
`estcell-all_individuals.EA-*-20260910` cells — record
`git_commit = f8fd01a74f2bcbf5f6824c52d4885a23e6db4eeb`, and that commit
contains **zero** files under `01b_estimation_cells/`:

```
$ git ls-tree -r f8fd01a74 --name-only | grep 01b_estimation_cells
$          # no output
```

`f8fd01a74` ("Fix stage 11 calling write_atomic outside the shared-helper env.",
2026-09-10) was HEAD while this module's `_h/` code still sat uncommitted in the
working tree, which is exactly what `git_dirty = true` is recording. It is not an
ancestor of the commit that landed the code, so the recorded commit is not merely
earlier than the code — it is on a different line of history.

**Use `de1795f90` instead** ("Add donor-group estimation cells: pooled discovery,
within-group SNP modeling.", 2026-09-10): the first HEAD-reachable commit
containing `01b_estimation_cells/_h/`, committed the same day the six runs were
produced.

That substitution is not a guess. The six runs' four recorded config checksums
reproduce **exactly** from `de1795f90` and not from the commit they record:

| field | recorded by the six runs | `config/` at `de1795f90` | `config/` at `f8fd01a74` |
|---|---|---|---|
| `config_cohorts_sha256` | `0021abb74695…` | `0021abb74695…` ✓ | `f2ecbb5c2291…` ✗ |
| `config_paths_sha256` | `02494c77256d…` | `02494c77256d…` ✓ | — |
| `config_thresholds_sha256` | `67ad135f3144…` | `67ad135f3144…` ✓ | — |
| `config_covariates_sha256` | `e5fe755de033…` | `e5fe755de033…` ✓ | — |

So the working tree at run time held `de1795f90`'s configuration while HEAD was
`f8fd01a74` — this module's code and its `cohorts.yml` estimation-cell block were
both uncommitted, which is what `git_dirty = true` records. `de1795f90` is the
commit that subsequently captured that tree.

It remains the closest recoverable code state rather than a proven-identical one:
`git_dirty = true` means the tree could still have differed from `de1795f90` in
files whose checksums are not recorded, and `_h/00_new_run.R` and
`_h/05_report_acceptance.R` were both substantially rewritten afterwards (at
`852ba936`). A regeneration from `de1795f90` should be checked against the sealed
outputs, not assumed to match them.

The three `estcell-AA.n118r{1,2,3}-caudate-20260918` cells record
`852ba93633918ce0a8ee172e99a582a310463ad8`, which **does** contain the module's
code. Only the six 2026-09-10 cells are affected.

**Why this is not a data defect.** Beyond the config checksums above, each cell's
loci verify as true subsets of the upstream Module 01 run by BED SHA-256 — a cell
copies loci verbatim and re-partitions donors over a fixed locus set, so its
scientific content is anchored to the upstream run rather than to this module's
code version. No VMR, phenotype, genotype or score is re-derived here.

**Recurrence is now detected at run time.** `00_shared/runid.R` records
`git_commit_has_module_code` in every new manifest and emits a loud
`PROVENANCE WARNING` when the recorded commit contains no file under the module
root. It warns rather than stops: whether a `git_dirty` run, or a run whose
commit lacks its own module's code, should be **refused outright** is a policy
decision for the PI (AGENTS.md §12), not an agent's to make. Until that is
decided, check the field before accepting a run.

## Regenerating a run

1. Check out the `git_commit` recorded in its `manifest.tsv` — but first read the
   section above, and check `git_commit_has_module_code` in the manifest. For the
   six 2026-09-10 cells that field is absent (it postdates them) and the recorded
   commit is unusable; use `de1795f90`.
2. Confirm the `config_*_sha256` fields match the current `config/`.
3. Re-run `_h/submit_estimation_cell.sh` from `01b_estimation_cells/_m`.

A completed SLURM job is not proof of scientific validity. Check
`task_reconciliation.tsv`, then `_h/05_report_acceptance.R`, then the module
README's acceptance gate.
