# 05_cpg_meqtl_burden — orthogonal genetic evidence

Asks whether increasing calibrated VMR local genetic control is associated with a greater fraction of constituent CpGs having conventional cis-meQTL support.

**Status: not implemented.** Gated on `01_vmr_catalog` and `02_local_genetic_variance` acceptance (AGENTS.md §6: "No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID").

## Migrating from

`meqtl-validation/01_cpg_meqtl_mapping/` and `meqtl-validation/02_vmr_meqtl_burden/`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Requirements

Use the corrected CpG-to-VMR membership keyed by the accepted `vmr_set_id`.
Report tested CpGs separately from prepared-but-untested CpGs. Model continuous
calibrated local variance as the primary predictor, with overdispersion-appropriate
and donor-robust inference. Matched extreme-group analysis is secondary evidence
only. Audit every concordance denominator, and resolve genomic inflation before
the figure freeze.

Internal meQTL mapping is convergent evidence, not independent replication.
Positive-only public resources cannot provide an external gradient, because
absence from a positive list is not a tested negative.

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
