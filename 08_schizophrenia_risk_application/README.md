# 08_schizophrenia_risk_application — required disease application

Tests whether schizophrenia-risk variants regulate methylation within genetically anchored VMRs. Intended for the main text, conditional on surviving corrected VMRs and the new local-genetic-control axis.

**Status: not implemented.** Gated on `05_cpg_meqtl_burden` and `06_transcription_splicing_coupling` acceptance (AGENTS.md §6: "No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID").

## Migrating from

`meqtl-validation/08_schizophrenia_risk_application/`, whose Phase 7 decision file records `retain_main_text_proof_of_application`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## That decision does not carry forward

The existing Phase 7 result (31 caudate loci, 361 pairs, 38 VMRs, eight
TX-coupled VMRs) was conditioned on the legacy predictability metric and
pre-repair VMR sets. Per AGENTS.md §8 it is a hypothesis to retest. Hero loci
`rs8048039` and `rs13331198` remain **candidates**, subject to reprioritization
after the corrected run.

## Preserved design principles

PGC schizophrenia loci are defined independently of methylation results.
Risk-variant-CpG tests keep their own FDR family. Association is tested against
the relative local SNP contribution score, not absolute PVE or legacy
predictability. External GTEx eQTL evidence
is support, not proof of mediation. At most five illustrative loci are
prioritized by a prespecified rule.

Add the integration analysis asking whether schizophrenia-linked VMRs are
enriched for LINE/L1, H3K9me3, quiescent chromatin, high-mappability repeat
intervals, and expression/splicing coupling. Positive connects the disease
application to the repeat/repressive architecture; null presents Phase 7 as a
separate proof of disease relevance.

Retention criteria and forbidden claims are in `config/schizophrenia.yml`.

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
