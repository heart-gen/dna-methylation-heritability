# 04_repeat_repressive_architecture — primary biological analysis

Tests whether a higher relative local SNP contribution score (`local_snp_contribution_score_z`, Module 02) is associated with repeat-rich and repressive genomic compartments. This is the module the manuscript's central claim rests on.

**Status: implemented, smoke-verified, not yet run in production.**
Stages 01-03 exist in `_h/` and have been run on a 400-VMR smoke in all three AA
regions (`rra-smoke-AA-{caudate,dlpfc,hippocampus}-20260823`). All 14 declared
covariates build at ~0 missingness, the hg38 -> hg19 liftover reports
11,223/11,341 uniquely mapped VMRs, and `03_apply_gates.R` emits
`permitted_claim` for all three outcomes. At 400 VMRs every outcome returns
"not supported: no region survives the locked sensitivities" -- that is the
expected result of a powerless smoke, **not** a scientific finding.

`config/repeat_annotations.yml` is `pi_locked: true`: Roadmap consolidated
epigenomes E068 (caudate), E073 (DLPFC), E071 (hippocampus), H3K9me3
gappedPeak plus the 15-state ChromHMM `15_Quies` state, both hg19, with hg38
VMRs lifted down using `inputs/supportfiles/_m/hg38ToHg19.over.chain`.

The chain is 01 features (per cell) -> 02 association (per cell) -> 03 gates
(cross-region) -> 04 figures + 05 seal (cross-region). Stages 03-05 span all
three regions by construction: the interpretation gates ARE the cross-region
rule, and no cell is sealed before the table that says what it may claim exists.
The three smoke cells sealed as
`GATES_APPLIED_0_OF_3_OUTCOMES_SUPPORTED_SMOKE_ONLY_NOT_ACCEPTABLE`.

This module consumes `03_local_snp_prediction` for its secondary predictor, so
its production driver is correctly refused by `require_accepted_upstream()`
until Module 03 records an accepted production run. That refusal was exercised:
`DRY_RUN=1 bash _h/run_repeat_architecture.sh AA` stops at
`00_new_run.R` with the AGENTS.md 6 message. The smoke runs above were created
with `--allow-unlocked`, which is the only path that bypasses it and which also
stamps `smoke_run = TRUE`.

## Migrating from

Repeat and cell-composition modules in `meqtl-validation/07_repeat_mappability_sensitivity/` and `meqtl-validation/11_celltype_compartment_sensitivity/`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Design

Primary outcomes: H3K9me3 overlap, quiescent-chromatin overlap, LINE/L1 overlap
or overlap fraction.

Primary predictor: standardized continuous
`local_snp_contribution_score_z`, constructed from the within-cell rank of the
frozen joint estimator among eligible loci. Secondary predictor: honest
`r2_pred_oof`, to ask whether the same
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

## Open work: the `exclude_snp_proximal` sensitivity is misimplemented

Opened 2026-09-02, after the first production run (`rra-AA-*-20260902`).

`config/repeat_annotations.yml` names this sensitivity
`sensitivities.exclude_snp_proximal_cpgs` -- a CpG-level exclusion, guarding
against the measurement artifact where a SNP under or beside a CpG destroys or
creates the site, so apparent methylation variance there is genotype rather than
epigenetics. `_h/02_test_association.R` realizes it at VMR level as
`d[snp_proximal_frac == 0]`: drop any VMR containing even one SNP-proximal CpG.

The median VMR has 71-74% of its CpGs SNP-proximal, so that cliff retains only
17.6-18.6% of loci, and the retained subset is biased (hippocampus):

| set        |    n | LINE/L1 | H3K9me3 | accessible | CpGs | length |
|------------|------|---------|---------|------------|------|--------|
| kept (arm) | 1634 |   0.047 |   0.050 |      0.429 | 13.5 | 583 bp |
| dropped    | 7631 |   0.069 |   0.071 |      0.366 | 15.3 | 768 bp |

The arm keeps short, CpG-poor, accessible, repeat-poor loci. Every outcome here
is an overlap FRACTION, and overlap with a LINE element or an H3K9me3 domain is
length-dependent, so the arm depletes the signal it is meant to test for reasons
unrelated to SNPs. Compounding this, common-SNP density is itself higher in
repeat-rich, late-replicating sequence, so conditioning on "no common SNP near
any constituent CpG" partly conditions on "not heterochromatin" -- removing the
exposure while testing for it.

Because `survives()` is a strict conjunction, this arm alone decides verdicts.
It blocked accessible-chromatin depletion in all three regions (primary
estimates -0.089 / -0.130 / -0.092, directional p to 3.7e-06, four of five arms
concordant) and dropped hippocampal LINE/L1 (primary +0.194, p=1.9e-04; this arm
+0.082, p=0.557), taking LINE/L1 from 2/2 to 1/2 regions.

### Planned, in order

1. **Graded VMR-level criterion**, as a quick review of the sensitivity: drop
   only loci DOMINATED by SNP-proximal CpGs (`snp_proximal_frac > 0.8`) rather
   than the `== 0` cliff. Targets artifact-driven VMRs, retains most loci, and
   largely removes the length bias. The current `== 0` arm is retained as a
   descriptive extreme case, reported but not gating.
2. **CpG-level recomputation**, the faithful realization: recompute each VMR's
   methylation summary from its non-SNP-proximal CpGs only, keeping the genomic
   span intact, then re-derive the predictor. The annotation outcome is
   unchanged, so no length bias is introduced. Requires per-CpG data from
   `01_vmr_catalog/_m/runs/{RUN}/cpg/chr_{N}/`, which this module does not
   currently read.

Amend as fidelity to `exclude_snp_proximal_cpgs`, NOT as relaxing a threshold.
The change favours the hypothesis, which is exactly when prespecification
matters; the defence is that the arm does not implement its own config key, and
that defence only holds if the amendment is written that way.

Neither step rescues the primary result: H3K9me3 and quiescent are null on the
UNRESTRICTED primary fits (|estimate| <= 0.087, p >= 0.13 in all three regions),
before any sensitivity applies. What is at stake is the LINE/L1 claim and the
accessible / H3K27ac depletion contrast.

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.

## Accepted runs

_(none)_

AGENTS.md §6 makes acceptance a human step: no row appears here until a
production run's gate stage passes and the PI records it. A smoke run is never
entered in this table.

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
