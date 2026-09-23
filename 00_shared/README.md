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
| `runid.R` | `new_run()`, `make_run_id()`, `seed_for()`, `write_atomic()`, `reconcile()`, `close_run()`, `seal_run_dir()` |
| `slurm.sh` | Repo-root resolution, conda env paths, `log_message()`, `run_r()`, `chrom_size()` |

## Why this exists

`vmr-analysis/` carried six near-identical copies of the same pipeline — two
cohort arms × three regions — with cohort and region as hard-coded string
literals. The copies drifted, and the drift *was* the defect register: V2 (region
filter), V3 (stale paths), V5 (function ignoring its argument), V6 (missing
bounds guard), V11 (incompatible cis windows), and the hippocampus copy reading
the DLPFC blacklist. "Do not duplicate region- or
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
- **Paths from config.** No Quest path inside an analysis function.
- **Unlocked config blocks production.** `assert_locked()` stops a run that would
  consume a PI decision nobody has made; smoke runs pass `--allow-unlocked` and
  get a loud warning.
- **Runs are immutable.** `new_run()` refuses to reuse a directory; `close_run()`
  makes it read-only, files *and* directories (see below).
- **Nothing vanishes quietly.** `reconcile()` requires every expected task to be
  completed, excluded, QC-failed, or failed, and refuses to close on an
  unexplained failure.

## Run sealing, and what it does not cover retroactively

`close_run()` calls `seal_run_dir()`, which clears the write bit on every file
**and every directory** in the run, preserving the setgid bit these trees carry.
Both halves are necessary: POSIX consults the *directory's* write bit when a
name is created or unlinked, never the file's own mode, so a run whose files are
all `0444` inside `2775` directories can still have any result deleted and
written back under the same name, silently and without disturbing
`output_checksums.tsv`.

Until 2026-09-23 only the files were chmodded, so that was true of every module
except `02_local_genetic_variance`, which never used `close_run()` and runs its
own `chmod -R a-w` in `_h/06_finalize_observed_run.R`. The asymmetry is how the
defect was found: in the run-retirement cleanup, Module 02 runs were the only
ones that resisted deletion.

**The repair binds new runs only.** Rewriting the modes of a run that has already
been closed and cited would itself be a modification after close, which AGENTS.md
5.2 forbids. Every run sealed before 2026-09-23 outside Module 02 therefore still
has writable directories. That is a property of those runs, not an outstanding
repair — treat their immutability as a convention rather than an enforcement, and
rely on `output_checksums.tsv` to detect a change rather than on the filesystem
to have prevented one.

## Tests

`tests/` is gitignored — hand-run before submitting an array:

```bash
Rscript -e 'source("00_shared/load.R"); testthat::test_dir("00_shared/tests")'
```

Covers: row-shuffle invariance, loud failure on missing and
duplicate donors, numeric chromosome ordering, the headerless-`.psam` case, run
immutability -- including that a sealed run genuinely refuses `create` and
`unlink`, keeps its setgid bit, and stops rather than reporting success when the
chmod does not take hold -- and refusal on unlocked config.
