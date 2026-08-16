# 06_transcription_splicing_coupling — regulatory consequences

Tests whether meQTL-supported or locally controlled VMRs are more likely to have existing significant associations with gene/transcript abundance or transcript usage/splicing.

**Status: not implemented.** Gated on `05_cpg_meqtl_burden` acceptance (AGENTS.md §6: "No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID").

## Migrating from

`meqtl-validation/06_transcription_splicing_integration/` and `meqtl-validation/09_libd_eqtl_mapping/`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Scope discipline

Do **not** initiate an unbounded transcriptome-wide fishing analysis. Reuse the
prespecified expression and splicing analyses and document their tested universe
explicitly. Adjust for the number of tested features, VMR length, VMR-to-feature
distance, methylation variance, local SNP number, and applicable technical
factors.

Allowed: "Genetically regulated VMRs are more frequently transcriptionally
coupled."
Forbidden: "Methylation mediates the genetic effect on expression or splicing."

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
