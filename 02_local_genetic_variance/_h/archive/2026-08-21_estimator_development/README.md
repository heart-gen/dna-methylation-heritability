# Archived Module 02 estimator-development scripts

These 61 files were moved out of the active `_h` root on 2026-08-21, after a
Quest `squeue` audit showed no active Module 02 job. No file was deleted.

They document the historical elastic-net calibrator, observed-EN runs,
high-signal audits, BSLMM comparisons, boundary recalibration, settings screen,
terminal joint-estimator experiment, and bounded recovery submissions. Their
scientific results remain documented in the module README and PI summary.

## Status

- `archived_not_production`: these paths must not be invoked for a new run;
- scripts retain their original filenames and version-control history;
- immutable `_m/runs/*` snapshots remain the authoritative executable record
  for completed historical runs;
- relative paths inside archived launchers are not maintained as an active API.

`ARCHIVE_MANIFEST.tsv` records the role and removal gate for every file.

## Later removal gate

Remove this directory only after all of the following are true:

1. the final joint training/validation and recovery runs are preserved in a
   tagged commit, archived release, or verified immutable run snapshot;
2. documentation no longer points to these files as executable entry points;
3. `MIGRATION_MANIFEST.tsv` records validated replacement by the active
   observed-score pipeline;
4. no downstream consumer or scheduler submission references an archived path.

Until then, this is a forensic archive—not an alternative pipeline.
