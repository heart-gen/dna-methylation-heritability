# 08_region_donor_generalization — cross-region robustness and identifiable regional heterogeneity

Establishes what **reproduces** across brain regions, which regional difference
is actually **identified**, and what the donor-group and matched-subset
contrasts can support.

**Status: not implemented.** Gated on `04_repeat_repressive_architecture` and `05_cpg_meqtl_burden` acceptance (AGENTS.md §6: "No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID").

## Migrating from

Donor-group, cross-region, and downsampling modules: `meqtl-validation/04_cross_region_sharing/`, `05_donor_group_comparison/`, `10_downsampling_caudate/`, `13_vmr_universe_nmatched/`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Scope, set by the PI 2026-09-06 before implementation

This module is **not** a symmetric "shared versus context-dependent" analysis.
Brain region is perfectly confounded with sequencing batch (AGENTS.md §8.1;
`writing-notes/WGBS_BATCH_REGION_CONFOUNDING.md` carries the evidence and the
full design review), so replication and difference are not two equal findings
here. Every output carries exactly one tier:

| tier | analysis | licenses |
|---|---|---|
| **Strongest** | cross-region replication | robustness of genetic-control architecture across technical and regional contexts |
| **Identified difference** | DLPFC vs hippocampus, both within batches 1–2 | genuine regional heterogeneity |
| **Mechanistic sensitivity** | caudate downsampled to n=118 | whether donor count explains the caudate excess — and nothing more |
| **Descriptive only** | caudate vs other regions | nothing; reported because readers will ask |

1. **Cross-region replication is the primary deliverable.** Describe it as
   robustness across technical and regional contexts. Do not write that the
   confounding strengthens the result: the batch boundary makes a successful
   replication more compelling, and does nothing for the interpretability of a
   difference. Existing material that belongs here and is currently only a
   by-product elsewhere: Module 03's check C, Module 05's 3% beta agreement
   across cells, Module 07's expression coupling.
2. **DLPFC vs hippocampus is the only clean region-difference analysis.** Use
   the Module 04 template: primary claim on the identified contrast,
   strict-conjunction sensitivity gating, confounded cell reported separately.
3. **Caudate downsampled to n=118** tests Module 03's untested attribution of
   the caudate excess to donor count (153 vs 118). Read it in two directions
   only — the excess largely disappears (donor count is a plausible major
   contributor) or it persists (donor count does not explain it). It **cannot**
   establish that a residual is biological; caudate stays batch-confounded
   whatever it shows. This limit belongs in the module's interpretation
   constraints, not in the reader's head.
4. **Caudate-vs-other differences stay descriptive.** Reuse the Module 04
   mechanism — `interpretation.technically_confounded_regions` sets the cell
   aside from the claim while keeping the estimate fitted, written and surfaced
   in dedicated columns.

### Donor-group axis

**PI decision 2026-09-06:** follow the v1 approach — compare AA and EA local
genetic-control estimates **on the common pooled-discovery VMR set**, the
`all_individuals` catalog, rather than on cohort-specific catalogs. Discovery
happens once in the pooled sample; the donor groups are then two disjoint sets
of donors evaluated on one fixed, shared locus set, so neither group's VMR
calling can advantage it.

- Do **not** contrast `AA` against `all_individuals`. Those are nested —
  `all_individuals` contains the AA donors — and a set-versus-superset
  comparison is not a donor-group contrast.
- **An EA estimation cell does not yet exist** in Modules 01 or 02. It must be
  defined and built before this axis can run.
- This axis is **not** exposed to the batch confounding: donor groups interleave
  within each region's libraries.
- Expect unequal n. Report it, and match or downsample as a tier-3 sensitivity
  rather than adjusting for it post hoc.

## Priorities

Prioritize biological generalization of local variance, repeat enrichment, meQTL
burden, and effect direction. Cross-population predictor portability is optional
and secondary.

AGENTS.md §7.6 already forbids raw score-level comparison across regions, so the
Module 02 score cannot be contrasted across cells at all. The region axis was
always narrower than the migration sources imply.

Do not attribute differences to ancestry-specific biology without eliminating
sample size, MAF, LD, SNP availability, assay, covariate, and brain-region
explanations. Use donor-group or population language approved by the PI:
`AA` = "Black American", `EA` = "non-Hispanic white American".

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
