# 06_partitioned_heritability — genome-wide disease relevance

Tests whether common-variant heritability for brain-relevant traits is enriched
in a **continuous** annotation built from the relative local SNP contribution
score (`local_snp_contribution_score_z`, Module 02), using stratified LD score
regression.

**Status: not implemented.** Gated on `02_local_genetic_variance` acceptance
(AGENTS.md §6: "No downstream production run may consume an upstream result
until the upstream README records a passing acceptance gate and immutable run
ID"). Module 02 has accepted runs for all six cells, so this module is
unblocked.

This module depends only on Module 02. Its position in the AGENTS.md §6
dependency list is a total order, not a claim that it consumes Modules 03–05.

## Why this module exists

`config/analysis_thresholds.yml` names `adds_nothing_beyond_nonsignificant_sldsc`
as an `omit_or_supplement_if` criterion for the schizophrenia application in
Module 09. Until this module produces a result, that criterion cannot be
evaluated. PI decision 2026-08-26: build S-LDSC as its own module rather than a
Module 09 sub-analysis, because it is SNP-level and genome-wide where every
other module is VMR-level.

## Migrating from

`local-snp-prediction/BA_only/tissue_comparison/clinical_enrichment/s-ldsc/`
(and the `all_individuals` counterpart), including `make_annot_continuous.py`,
`ldsc_wrapper.py`, `region_heritability.py`, `fdr_correction.py`, and
`interpreting_sldsc_results.md`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Endpoint discipline

The annotation is **continuous**. Do not rebuild the retired v1 form, which
partitioned VMRs into heritable and non-heritable classes at a threshold on
`r_squared_cv` — banned by AGENTS.md §3. No threshold, no grouping, and no
absolute locus PVE enters the annotation.

Report the metrics named in the legacy `interpreting_sldsc_results.md`:
enrichment, enrichment p, and the tau coefficient z-score. A significant
enrichment on a continuous annotation is a statement about where common-variant
heritability concentrates, not about any individual VMR.

## Unresolved before production

**LD reference for an AA primary arm.** Standard S-LDSC LD scores are
EUR-derived, while the primary cohort is Black American donors and the
annotation is defined on VMRs discovered in that cohort. v1 carried separate
`_AA` and `_EA` variants of `region_heritability.py`. The defensible position is
that the annotation is a genomic feature rather than a donor property, but that
must be stated explicitly in Methods rather than assumed. Resolve before the
first production run.

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
