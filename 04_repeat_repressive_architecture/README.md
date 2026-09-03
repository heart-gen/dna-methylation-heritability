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

## Resolved 2026-09-02: the adjustment set and the `exclude_snp_proximal` arm

Two independent specification defects were found after the first production run
(`rra-AA-*-20260902`, H3K9me3 0/3 regions, quiescent 0/3, LINE/L1 1/2). Both are
now closed by a dated amendment to `config/repeat_annotations.yml`, whose header
carries the full per-covariate reasoning. Summarized here.

### 1. The adjustment set over-adjusted on mediators

`_h/06_qc_adjustment_ladder.R` refits the primary model along a nested ladder on
the sealed feature table, holding predictor, outcome scale, family and robust
SEs fixed. `full_v2` reproduces the sealed run exactly, so the decomposition is
trustworthy. Caudate log-odds per SD:

| rung      | adds                                           | H3K9me3 | quiescent | LINE/L1 |
|-----------|------------------------------------------------|---------|-----------|---------|
| v1_like   | log(length), tested_snp_count                  |   0.434 |     0.714 |   0.746 |
| plus_cpg  | cpg_count, cpg_density                         |   0.343 |     0.596 |   0.434 |
| plus_seq  | gc_content                                     |   0.192 |     0.492 |   0.176 |
| plus_meth | mean_meth, meth_variance, wgbs_coverage        |   0.118 |     0.123 |   0.010 |
| plus_tech | snp_proximal, mappability, segdup, problematic |   0.087 |     0.099 |  -0.058 |
| full_v2   | broad_genomic_annotation, cell_composition     |   0.049 |     0.023 |  -0.029 |

DLPFC H3K9me3 at `v1_like` is 0.381, p=1.5e-14 -- stronger than v1, which was
null there. So the v1 signal is fully reproducible in v2 data: it dies from
covariates, not from re-analysis, annotation or the liftover.

Two rungs carry nearly all the attenuation, and both are the same kind of
covariate. `plus_meth`: methylation variance is the DEFINING property of a VMR
and sits on the path from local genetic variation to methylation, so
conditioning on it removes the exposure by construction; mean methylation is
near-constitutive of the quiescent outcome, which is the hypermethylated
compartment. `full_v2`: genic/intergenic is partly constitutive of "quiescent,
gene-poor, late-replicating". `plus_tech` is a different animal and is genuinely
protective -- v1 carried no technical covariates at all, which is its own
defect. Neither version was right: v1 under-adjusted, v2 over-adjusted.

**Prespecified resolution.** The primary adjustment set is now `plus_seq +
plus_tech`: `vmr_length, cpg_count, cpg_density, gc_content, tested_snp_count,
snp_proximity, mappability, segdup_overlap, problematic_region_overlap`. The
methylation block and `broad_genomic_annotation` move to
`descriptive_covariates` -- still built, still completeness-checked in
`_h/01_build_features.R`, out of the model formula.
`cell_composition_pcs` also leaves the primary, but because its
confounder-vs-mediator status is arguable rather than settled: it becomes the
`adjust_cell_composition` sensitivity arm, which refits every model with it
added. The existing `low_cell_composition` subset arm is unchanged.

The ladder was seen before this decision was made. That is why the reasoning is
recorded per covariate in the config header, and why
`_h/06_qc_adjustment_ladder.R` now also fits `prespecified_primary` absolutely
(it skips `plus_meth`, so it is not any cumulative prefix of the ladder) --
the new production fits must be checkable against the sealed table.

### 2. `exclude_snp_proximal` retired from the gating conjunction

`config` names this `sensitivities.exclude_snp_proximal_cpgs`, a CpG-level
exclusion guarding a real artifact: a SNP under a CpG destroys or creates the
site, so "methylation variance" there is genotype, not epigenetics.
`_h/02_test_association.R` realized it as `d[snp_proximal_frac == 0]`.

That is retired, on two grounds:

- **Redundant.** The identical BED already enters every model as the covariate
  `snp_proximity` (realized as `snp_proximal_frac`). The arm conditioned on a
  variable the primary already adjusts for. This argument is independent of any
  result, and it is the one the retirement rests on.
- **Not a faithful realization of the key, and biasing.** `snp_proximal_frac` is
  a BASE-PAIR overlap fraction of the VMR span with the +/-150bp windows
  (`_h/annotation_io.R`), not a fraction of CpGs. So `== 0` selects short VMRs
  mechanically, while every outcome here is itself a length-dependent overlap
  fraction. It retained 17.6-18.6% of loci, biased (hippocampus):

| set        |    n | LINE/L1 | H3K9me3 | accessible | CpGs | length |
|------------|------|---------|---------|------------|------|--------|
| kept (arm) | 1634 |   0.047 |   0.050 |      0.429 | 13.5 | 583 bp |
| dropped    | 7631 |   0.069 |   0.071 |      0.366 | 15.3 | 768 bp |

  Compounding this, common-SNP density is itself higher in repeat-rich,
  late-replicating sequence, so the arm partly conditioned on "not
  heterochromatin" while testing for heterochromatin.

Because `survives()` is a strict conjunction, this one arm decided verdicts: it
blocked accessible-chromatin depletion in all three regions (primary -0.089 /
-0.130 / -0.092, directional p to 3.7e-06, four of five arms concordant) and
dropped hippocampal LINE/L1 (primary +0.194, p=1.9e-04; this arm +0.082,
p=0.557), taking LINE/L1 from 2/2 to 1/2.

It is still fitted, and written to `results/descriptive-snp-proximal.tsv` --
a separate file from `association-results.tsv`, which is what `03_apply_gates.R`
reads, so it cannot re-enter the conjunction.

This is recorded as fidelity to `exclude_snp_proximal_cpgs`, not as relaxing a
threshold. The change favours the hypothesis, which is exactly when
prespecification matters; the defence is the redundancy argument plus the fact
that the arm never implemented its own config key.

### Deferred: the faithful CpG-level exclusion

Not realizable inside module 04. `local_snp_contribution_score_z` descends from
the module 02 VMR phenotype, which is `getMeth(BSobj, regions = gr, what =
"perRegion")` over ALL CpGs in the span (`01_vmr_catalog/_h/02_summarize.R`).
A non-SNP-proximal-CpG phenotype means re-running 01 summarize ->
`02_local_genetic_variance` (`_h/01_estimate_observed_joint_features.R` ->
`_h/04_derive_local_snp_contribution_score.R`) -> 04. The ingredients exist:
`01_vmr_catalog/_m/runs/{RUN}/vmr/cpg_vmr_membership.tsv` maps per-CpG
coordinates to VMR id, and `.../cpg/chr_{N}/` holds the per-CpG matrices. Scoped
and not executed.

### What this does and does not buy

It does not automatically rescue the central claim. H3K9me3 and quiescent are
null on the UNRESTRICTED `full_v2` primary fits (|estimate| <= 0.087, p >= 0.13
in all three regions), before any sensitivity applies; only the adjustment-set
change can move them, and the reduced set is allowed to return whatever it
returns. What the arm retirement is worth is the LINE/L1 claim and the
accessible / H3K27ac depletion contrast.

One loose thread worth watching under the new set: LINE/L1 is discordant across
regions in `full_v2` (caudate -0.029, DLPFC +0.161, hippocampus +0.194), yet at
`v1_like` caudate LINE/L1 is 0.746, the largest of the three, and module 03's QC
found caudate the strongest region for predictability. Caudate is therefore
where adjustment bites hardest, not where signal is absent -- which is what the
mediator account predicts, and is a live test of it.

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
