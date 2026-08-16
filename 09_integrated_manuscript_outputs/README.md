# 09_integrated_manuscript_outputs — one source of truth

Consumes only accepted immutable upstream runs and produces every manuscript number, table, and figure.

**Status: not implemented.** Gated on acceptance of every upstream module it consumes (AGENTS.md §6: "No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID").

## Migrating from

Manuscript figures and consolidated tables currently scattered across `meqtl-validation/12_supplementary_data/` and per-module figure directories.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Products

Manuscript-number registry; main and supplementary tables; main and supplementary
figures; figure source-data tables; analysis-to-claim matrix; exclusions and
denominator table; software and run manifest; manuscript-ready Methods and
Results summaries.

## Rule

Do not manually assemble final figures from files copied across old directories.
Every figure panel must record its source run ID, table, script, and filter.

Always report denominators, exclusions, brain region, donor group, VMR set, and
the exact metric used.

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
