# 04_repeat_repressive_architecture — primary biological analysis

Tests whether increasing calibrated local genetic control is associated with repeat-rich and repressive genomic compartments. This is the module the manuscript's central claim rests on.

**Status: not implemented.** Gated on `02_local_genetic_variance` acceptance (and `03_local_snp_prediction` for the secondary predictor) (AGENTS.md §6: "No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID").

## Migrating from

Repeat and cell-composition modules in `meqtl-validation/07_repeat_mappability_sensitivity/` and `meqtl-validation/11_celltype_compartment_sensitivity/`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Design

Primary outcomes: H3K9me3 overlap, quiescent-chromatin overlap, LINE/L1 overlap
or overlap fraction.

Primary predictor: standardized continuous `h2_en_calibrated` among interpretable
loci. Secondary predictor: honest `r2_pred_oof`, to ask whether the same
compartments are enriched among *imputable* VMRs.

Adjustment/matching variables and locked sensitivities are enumerated in
`config/repeat_annotations.yml`.

## Interpretation gates

- H3K9me3 and quiescent enrichment may be described as shared across regions only
  if direction and inference survive locked sensitivities in **all three**.
- LINE/L1 may be described as multi-region if it survives in **at least two**. If
  it survives only in caudate, it is presented as caudate-specific.
- A DLPFC reversal or null after the high-mappability restriction is never hidden.

Overlap alone never implies activity, expression, or retrotransposition.

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
