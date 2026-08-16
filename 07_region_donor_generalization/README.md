# 07_region_donor_generalization — boundaries of the biology

Establishes what is shared and what is context-dependent across brain regions, donor groups, and matched subsets.

**Status: not implemented.** Gated on `04_repeat_repressive_architecture` and `05_cpg_meqtl_burden` acceptance (AGENTS.md §6: "No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID").

## Migrating from

Donor-group, cross-region, and downsampling modules: `meqtl-validation/04_cross_region_sharing/`, `05_donor_group_comparison/`, `10_downsampling_caudate/`, `13_vmr_universe_nmatched/`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Priorities

Prioritize biological generalization of local variance, repeat enrichment, meQTL
burden, and effect direction. Cross-population predictor portability is optional
and secondary.

Do not attribute differences to ancestry-specific biology without eliminating
sample size, MAF, LD, SNP availability, assay, covariate, and brain-region
explanations. Use donor-group or population language approved by the PI:
`AA` = "Black American", `EA` = "non-Hispanic white American".

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
