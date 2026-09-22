# supplementary_data — the deposition list

One file is the point of this directory: `supplementary_data_manifest.tsv`. It
records, for every Supplementary Data item the manuscript will cite, **where the
real file already lives** in this repository.

Nothing is copied here and nothing is symlinked here. A staging tree would be a
second copy of ~420 MB that goes stale the first time a run is superseded, and a
tree of symlinks breaks the moment it is archived or moved off Quest. The
manifest points at the immutable run outputs in place, so the list can be
verified against the pipeline at any time instead of being trusted.

```
python3 supplementary_data/_h/verify_manifest.py           # report
python3 supplementary_data/_h/verify_manifest.py --write   # refresh sizes
```

The checker expands each location, counts the files, sums the bytes, asks git
whether they are tracked, and exits non-zero if a location that claims to be
deposit-ready resolves to nothing. Run it before staging an upload.

## Columns

| column | meaning |
|---|---|
| `sd` | Supplementary Data number, as the manuscript cites it |
| `title` | the item's title; repeated on every row of the same `sd` |
| `module` | the module that produced it |
| `run_ids` | the **accepted** run IDs, in the brace shorthand of the module READMEs |
| `location` | repo-relative path, glob, or `dir/{a,b}.tsv` brace group |
| `n_files`, `bytes` | computed by the checker, never edited by hand |
| `in_git` | `yes` / `partial` / `no`; see the tiering below |
| `distribution` | `git+zenodo`, `zenodo`, `quest`, `none` |
| `status` | see the table below |
| `donor_ids` | whether the file carries a donor or sample identifier |
| `notes` | what the file holds and what may not be claimed from it |

## Tiering

`_m/runs/` is gitignored, so a run output always answers `no` to `in_git`; only
`_m/combined/` is tracked. That gives three tiers:

- **`git+zenodo`** — a cross-region summary in `_m/combined/`. Readable from a
  clone, and mirrored to Zenodo so the DOI is self-contained.
- **`zenodo`** — a per-VMR or per-locus table from an accepted run. Too large for
  git, small enough to deposit, and generally the sufficient statistic for a
  reported number: without it a reader cannot recompute an FDR or refit an axis
  model, because the significant subset does not carry the family it was
  corrected within.
- **`quest`** — bulk intermediates that no claim rests on: full nominal meQTL
  output, per-fold diagnostics, sharded copies of a table deposited whole.

## Statuses

| status | meaning |
|---|---|
| `ready` | accepted run, no blocker; upload as-is |
| `pending_acceptance` | run is sealed and passes its gate, but has no acceptance row |
| `pending_build` | the output does not exist yet |
| `needs_deid` | depositable only after a donor or sample identifier is removed |
| `partial_deposit` | only a named subset ships; the resolved size overstates it |
| `transform_on_deposit` | ships, but a column must be rewritten when the archive is cut |
| `audit_only` | retained inside its module for audit; not a reportable quantity |
| `excluded` | deliberately Quest-only |
| `superseded` | a v1 output a v2 run replaces; retire per AGENTS.md §3 |
| `external` | not generated here; cite the source instead of redepositing |
| `blocked_donor_ids` | carries donor identifiers **and is tracked in git** |
| `blocked_unlocked_config` | produced by a stage whose config is not PI-locked |

Only `ready`, `needs_deid`, `pending_acceptance` and `transform_on_deposit`
count toward the deposit total the checker prints.

## Identifiers

`BrNum` is a de-identified donor identifier and may be retained in deposited
tables (PI decision, 2026-09-20). The `donor_ids` column stays in the manifest as
a factual record of which files carry one, not as a blocker. Two notes:

- The array barcode **is** stripped. `00_shared/identity.R::strip_sample_barcode`
  drops the `::barcode` suffix, and it refuses rather than proceeding if the strip
  would collapse two distinct donors into one row. It is applied in
  `03_combine_oof.R` (so future runs never emit it) and in
  `07_stack_donor_group_predictions.R`, after that stage's disjointness check has
  run on the full identifier. The three `_m/combined/` tables were regenerated on
  2026-09-20 and are byte-identical to the previous version with the suffix
  removed. The three `lsp-AA-*-20260825` runs were sealed before the change and
  are immutable, so their copies still carry it and are marked
  `transform_on_deposit`: strip at archive time, never by editing a sealed run.
- `01_vmr_catalog/.../vmr/phenotypes/` stays excluded regardless. Those are
  individual-level methylation measurements, not an identifier question, and
  access is governed by dbGaP phs000979.

## Still open before submission

**Module 10 is sealed but unaccepted.** SD13's tables all carry
`citable = FALSE`. The Zenodo record is versioned and public, so nothing from
`env-AA-*-20260920-a` should go on it before the acceptance rows exist.

**Module 11 has no accepted run.** SD14 is `pending_build`, and Figures 1-2 were
last built on run IDs retired on 2026-09-17.

**Resolved 2026-09-20:** the retired v1 elastic-net calibration outputs
(`02_local_genetic_variance/_m/combined/calibrated-local-h2-all-cells.tsv` and
`observed-run-qc-all-cells.tsv`) are untracked from git and gitignored under
AGENTS.md §3. They remain on Quest for audit and are recoverable from history.

## Relationship to the manuscript repository

`content/sdata/supplementary_data_manifest.tsv` in the manuscript repository is
the **distribution** ledger: archive name, bytes, SHA-256, and where each item is
served from. This file is the **provenance** ledger: which accepted run produced
each item and which file on disk it is. The `sd` numbers are the join key, and
they must be changed in both places at once.
