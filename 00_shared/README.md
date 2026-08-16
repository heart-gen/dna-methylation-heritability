# 00_shared — v2 library code

Shared utilities for the numbered v2 revision modules. This directory has **no**
`_h/` and no `_m/`: it produces no generated output, so the code/output split
that `_h/`/`_m/` exists to enforce does not apply. Files sit at the top level.

Load it at the top of any v2 script:

```r
source(file.path(Sys.getenv("V2_REPO_ROOT", "."), "00_shared", "load.R"))
```

```bash
source "$REPO_DIR/00_shared/slurm.sh"
```

## Files

| File | Purpose |
|---|---|
| `load.R` | Bootstrap; self-locates the repo root and sources the rest in order |
| `config.R` | `load_config()`, `resolve_path()`, `assert_locked()`, `parse_v2_args()`, `cohort_def()`, `sample_blacklist()` |
| `identity.R` | `align_by_id()`, `read_psam()`, `assert_no_dups()`, `assert_expected_n()`, `donor_checksum()` |
| `chrom.R` | `chrom_order()`, `sort_chunks()`, `verify_fid_iid()`, `sort_genomic()`, `has_ct_mask()` |
| `runid.R` | `new_run()`, `make_run_id()`, `seed_for()`, `write_atomic()`, `reconcile()`, `close_run()` |
| `slurm.sh` | Repo-root resolution, conda env paths, `log_message()`, `run_r()`, `chrom_size()` |

## Why this exists

`vmr-analysis/` carried six near-identical copies of the same pipeline — two
cohort arms × three regions — with cohort and region as hard-coded string
literals. The copies drifted, and the drift *was* the defect register: V2 (region
filter), V3 (stale paths), V5 (function ignoring its argument), V6 (missing
bounds guard), V11 (incompatible cis windows), and the hippocampus copy reading
the DLPFC blacklist. AGENTS.md §5.3: "Do not duplicate region- or
cohort-specific copies of the same code."

The single most important function here is `align_by_id()`. It exists because of
defect **V1**: the legacy code reordered a response matrix with `match()` while
merely subsetting the design matrix with `%in%`, so every donor's methylation
was regressed against a different donor's principal components. The shapes were
right, so nothing complained — for the whole life of the project. `align_by_id()`
reorders both sides and asserts the IDs match afterwards.

## Conventions it enforces

- **Identity before arithmetic.** Never subset two things and assume they line
  up. Missing or duplicated donors are errors, not dropped rows.
- **Paths from config.** No Quest path inside an analysis function
  (AGENTS.md §9).
- **Unlocked config blocks production.** `assert_locked()` stops a run that would
  consume a PI decision nobody has made; smoke runs pass `--allow-unlocked` and
  get a loud warning.
- **Runs are immutable.** `new_run()` refuses to reuse a directory; `close_run()`
  makes it read-only.
- **Nothing vanishes quietly.** `reconcile()` requires every expected task to be
  completed, excluded, QC-failed, or failed, and refuses to close on an
  unexplained failure.

## Tests

`tests/` is gitignored — hand-run before submitting an array:

```bash
Rscript -e 'source("00_shared/load.R"); testthat::test_dir("00_shared/tests")'
```

Covers AGENTS.md §10.1: row-shuffle invariance, loud failure on missing and
duplicate donors, numeric chromosome ordering, the headerless-`.psam` case, run
immutability, and refusal on unlocked config.
