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
- LINE/L1 may be described as multi-region if it survives in **at least two**
  ELIGIBLE regions. Since 2026-09-06 caudate is not eligible for this outcome --
  it is sequencing batch 3 and region is perfectly confounded with batch -- so
  the two eligible regions are DLPFC and hippocampus and both must survive. The
  caudate estimate is still fitted, still reported, and labelled technically
  uncertain. See "Resolved 2026-09-06: what LINE/L1 is allowed to claim".
- A DLPFC reversal or null after the high-mappability restriction is never hidden.
- The matched-measurability arm (`_h/08`) is secondary evidence and never gates.

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

### The caudate LINE/L1 collapse: attributed 2026-09-06

Under the prespecified set LINE/L1 is 0.038 (p=0.41) in caudate against 0.298
(DLPFC) and 0.330 (hippocampus), despite caudate having the LARGEST `v1_like`
estimate of the three (0.746). An earlier draft of this README read that as
"caudate is where adjustment bites hardest, not where signal is absent". That
reading is **wrong** and is retracted. `_h/07_qc_covariate_attribution.R`
attributes the collapse, on the sealed `rra-AA-*-20260902` tables:

- **`gc_content` is the single operative covariate, and it acts in opposite
  directions by region.** Added alone to `v1_like` it takes caudate 0.746 ->
  0.217, while it RAISES DLPFC (0.423 -> 0.496) and hippocampus (0.437 ->
  0.521). Dropped from the prespecified set it restores caudate to 0.259
  (p=4.1e-08) but lowers the other two. `wgbs_coverage` adds nothing beyond it
  in caudate (0.217 vs 0.219); it is a proxy, not a second cause.
- **The entanglement is in the PREDICTOR, not the annotation.** Spearman
  correlation of `local_snp_contribution_score_z` with `gc_content` is -0.477 in
  caudate against -0.035 and -0.012. LINE/L1 annotation is region-invariant, so
  this is a property of how the caudate score was constructed.
- **`cor(gc_content, wgbs_coverage)` is -0.807 in caudate and +0.529 / +0.515 in
  the other two** -- a sign flip on a basic assay/sequence relationship, not a
  magnitude difference. `wgbs_coverage` is read straight off the BSseq objects
  (`_h/01_build_features.R`), so this is not a Module 04 join artifact.

- **It is not specific to LINE/L1.** An earlier note here implied it was; that
  is also retracted. Every caudate outcome moves when `gc_content` is dropped
  (H3K9me3 0.123 -> 0.258, quiescent 0.442 -> 0.521, H3K27ac -0.422 -> -0.490).
  The bias on any outcome scales as `cor(score, gc) x cor(gc, outcome)`, and
  LINE/L1 simply has the largest product (0.224 in caudate, <=0.009 in the other
  two regions) because L1 is the most AT-rich annotation in the panel. LINE/L1
  is where the attenuation crosses zero, not where it uniquely occurs.

`gc_content` is in `plus_seq`, a retained and uncontested confounder. So the
caudate null is the adjustment set behaving correctly, and it is **not**
evidence for or against the mediator account, which concerns `plus_meth` and
`broad_genomic_annotation` only. What the caudate cell is allowed to say is
settled separately, in "Resolved 2026-09-06" below.

Upstream provenance is parallel across regions (`vmrcat-AA-*-20260816` ->
`lgv-AA-*-20260823` -> `lsp-AA-*-20260825`), so the asymmetry is in the caudate
data, not in the wiring. Caudate has 153 donors against 118 and 117, calls
11,530 VMRs against 9,572 and 9,497, and is sequenced shallower (mean per-CpG
depth 14.4 against 19.8). See "Open: the caudate coverage/GC inversion".

### Rejected 2026-09-06: a caudate-specific coverage filter or covariate set

Tested directly, not argued. The `coverage` pass of
`_h/07_qc_covariate_attribution.R` sweeps a minimum `wgbs_coverage` threshold;
if the caudate entanglement were a marginal-coverage artifact, a stricter filter
should relax it. It does the opposite:

| caudate `wgbs_coverage >=` | n | cor(score, gc) | cor(gc, cov) | LINE/L1 prespecified |
|---|---|---|---|---|
| 0 (as run) | 11341 | -0.477 | -0.807 |  0.038 (p=0.41) |
| 10 | 8957 | -0.509 | -0.787 |  0.034 (p=0.48) |
| 12 | 6085 | -0.539 | -0.750 | -0.003 (p=0.95) |
| 15 | 3442 | -0.504 | -0.676 | -0.007 (p=0.91) |
| 20 | 1515 | -0.312 | -0.369 |  0.002 (p=0.98) |

The score/GC entanglement gets STRONGER on better-covered loci, the coverage
inversion persists across the whole range, and the LINE/L1 estimate stays at
zero at every threshold. The drift toward zero correlation at `>= 20` is range
restriction on the thresholded variable at n=1515, not a recovery. DLPFC and
hippocampus meanwhile lose signal monotonically under the same filter (DLPFC
0.298 -> 0.116 at `>= 20`), so raising the Module 01 threshold from 5 to 10
would cost the two regions that carry the claim and buy caudate nothing.

The `set_overlap` pass rules out the VMR set as the source: the inversion is
present in caudate VMRs shared with both other regions (cor(gc, cov) = -0.769,
n=3416) and in caudate-unique ones (-0.802, n=7925) alike. It is a whole-region
property, so it is not the ~2,000 extra caudate VMRs.

Two conclusions follow. First, a stricter coverage filter is not the fix -- the
caveat on the screen is that it drops whole VMRs by mean depth rather than
re-calling them per CpG as `01_vmr_catalog/_h/00_prepare.R` would, so it is a
screen and it is the NEGATIVE result that is informative. Second, a
caudate-specific covariate set is not warranted: `gc_content` is already in the
prespecified set and already removes the entanglement, so there is no covariate
missing for caudate. Fitting a different adjustment set per region on the
strength of these numbers is the same after-the-fact selection the 2026-09-02
amendment was written to prevent, and it would make the three regions
non-comparable. The covariates stay locked and identical across regions.

### Resolved 2026-09-06: the inversion is a sequencing batch confounded with region

The standing hypothesis in an earlier draft of this section -- "a
caudate-specific library or conversion batch effect, untested" -- is now tested
and confirmed, and the earlier suggestion that a caudate re-export might help is
**retracted**.

These WGBS data are the AANRI phase 1 delivery
(<https://github.com/LieberInstitute/aanri_phase1>). Its `raw-data/batch-1` and
`raw-data/batch-2` carry 150 DLPFC and 150 hippocampus library tokens and **zero
caudate**. Caudate is `raw-data/batch-3`, and
`batch-3/bs_objs/batch3_combined/CpGassays.h5` (83,183,223,042 bytes) is
byte-size-identical to the `wgbs/caudate/raw/CpGassays.h5` that the caudate
BSobj build reads.

**Brain region and sequencing batch are perfectly confounded in this delivery.**
No contrast available in these data can separate them.

Evidence that the inversion is a library property and not an analysis artifact:

- Genome-wide on chr22 -- 1 kb bins, all CpGs, reference GC, no VMRs involved --
  `cor(gc, coverage)` is -0.638 to -0.719 in caudate against +0.678 (DLPFC) and
  +0.655 (hippocampus).
- All 308 caudate samples are individually negative, uniform across all seven
  brain-bank `source` values. There is no within-caudate batch.
- Donor identity is not the problem: `01_vmr_catalog/_h/00_prepare.R` resolves
  donors through `colData(BSobj)$brnum` with `assert_no_dups`, never by column
  position, despite caudate's `BSobj` having NULL `colnames`.
- Re-processing cannot fix it. `BSmooth`, blacklist removal, batch combination
  and BrNum mapping all leave the `Cov` assay untouched, so a GC-coverage curve
  set by PCR amplification and post-bisulfite size selection is not recoverable
  downstream.

This reaches past Module 04, so it is written up repo-level in
[`writing-notes/WGBS_BATCH_REGION_CONFOUNDING.md`](../writing-notes/WGBS_BATCH_REGION_CONFOUNDING.md),
which carries the per-module exposure assessment and the rules for what a
cross-region claim may say. Read it before writing any caudate-vs-other-region
sentence.

### Resolved 2026-09-06: what LINE/L1 is allowed to claim

Recorded as a dated amendment to the PI-locked `config/repeat_annotations.yml`.
The decision keeps every region in every analysis and changes only the evidence
base for one outcome:

| | decision |
|---|---|
| Main methylation / genetic-control analyses | keep all three regions |
| Region as a biological variable | unchanged; brain-region biology is real |
| `gc_content` adjustment | **keep**, but it is not treated as proof the LINE/L1 bias is solved |
| LINE/L1 inference | estimated region-specifically, as before |
| Primary LINE/L1 biological claim | rests on DLPFC and hippocampus, the two regions with no GC-coverage inversion |
| Caudate LINE/L1 | reported separately and labelled **technically uncertain**; neither corroboration nor refutation |

Realized as `interpretation.technically_confounded_regions: {line_l1_frac:
[caudate]}`. `03_apply_gates.R` counts survival over ELIGIBLE regions only, so
`line_l1_multiregion_requires: 2` now means both of the two eligible regions,
and the caudate estimate and p travel with the claim in the
`regions_excluded` / `excluded_estimate` / `excluded_p` columns of
`interpretation-claims.tsv`. It cannot be dropped by summarizing the table.

This is a **restriction** of the evidence base, not a relaxation: it removes a
region from a claim rather than admitting one.

**Also rejected 2026-09-06: dropping `gc_content`.** The proposal was that GC
matters only marginally in DLPFC and hippocampus while being strongly entangled
with coverage in caudate. The first half does not hold for LINE/L1, which is the
outcome at issue:

| region | LINE/L1 with `gc_content` | without |
|---|---|---|
| caudate | 0.038 (p=0.41) | 0.259 (p=4.1e-08) |
| DLPFC | 0.298 (p=6.4e-11) | 0.190 (p=1.7e-05) |
| hippocampus | 0.330 (p=4.7e-14) | 0.206 (p=8.6e-07) |

Dropping it costs DLPFC 36% and hippocampus 38% of a real effect while
manufacturing caudate significance out of the batch confound. It is marginal for
the other outcomes in those two regions, but not for LINE/L1 -- L1 is the most
AT-rich annotation in the panel, so it is exactly where GC does the most work.
And in the two clean regions GC adjustment RAISES the estimate; genuine
over-adjustment would lower it, which is evidence against reading GC as a
mediator here. `gc_content` is fixed reference sequence, upstream of both the
exposure and the outcome. It stays.

### Added 2026-09-06: the matched-measurability sensitivity

`_h/08_matched_measurability.R`, config
`sensitivities.matched_measurability` (which realizes the long-declared,
previously unimplemented `matched_high_low_comparison` key). **Secondary
evidence, non-gating**, written to its own files so `survives()` never sees it.

Regression adjustment for GC and coverage assumes a functional form and
extrapolates across the whole covariate range -- including where one group has
almost no support, which is precisely where a batch artifact lives. Matching
asks the cleaner question instead: are these loci different from loci that were
**equally measurable**? Two transposes, both reported:

- `annotation_vs_matched` -- annotation-overlapping loci vs non-overlapping loci
  matched on `gc_content`, `wgbs_coverage` and `vmr_length`; contrast on the
  predictor. The literal form of the question.
- `score_high_vs_low` -- top-tertile vs matched bottom-tertile predictor loci;
  contrast on the annotation fraction. The transpose that lines up with the
  primary model.

Greedy 1:1 nearest-neighbour on Mahalanobis distance of the standardized match
variables, no replacement, 0.2 pooled-SD caliper on *every* variable, seed and
parameters fixed in config. `MatchIt`/`optmatch` are not installed in the
project environment, so this is implemented directly and the greedy order is
seeded rather than left to row order. Balance is written beside the estimates in
`descriptive-matched-balance.tsv`; a matched analysis without a balance table is
not interpretable.

Run against the sealed `rra-AA-*-20260902` feature tables (LINE/L1):

| region | contrast | pairs | cases unmatched | estimate | p | max post-match \|SMD\| |
|---|---|---|---|---|---|---|
| caudate | annotation_vs_matched | 1012 | 649 / 1661 (39%) | 0.130 | 7.6e-05 | 0.013 |
| DLPFC | annotation_vs_matched | 1200 | 136 / 1336 (10%) | 0.429 | 3.0e-26 | 0.008 |
| hippocampus | annotation_vs_matched | 1253 | 123 / 1376 (9%) | 0.401 | 1.6e-24 | 0.005 |
| caudate | score_high_vs_low | 1773 | 2008 / 3781 (53%) | 0.018 | 7.0e-04 | 0.075 |
| DLPFC | score_high_vs_low | 1819 | 1295 / 3114 (42%) | 0.070 | 1.1e-22 | 0.032 |
| hippocampus | score_high_vs_low | 1830 | 1259 / 3089 (41%) | 0.078 | 1.5e-28 | 0.037 |

Three things to read off it, stated plainly because two of them cut against the
regression result:

1. **DLPFC and hippocampus survive the design change.** The primary LINE/L1
   claim does not depend on the functional form of the GC adjustment.
2. **The pre-match imbalance quantifies the batch.** Caudate L1-overlapping loci
   differ from non-overlapping loci by SMD -1.77 on GC and **+1.21** on
   coverage; DLPFC and hippocampus are -0.84 / -0.37 and -0.81 / -0.33. The
   coverage imbalance is opposite in sign, and caudate loses 39-53% of its cases
   to the caliper against 9-42% elsewhere. That is direct evidence that the
   caudate regression was extrapolating off common support.
3. **On common support, caudate is positive, not null** -- 0.130 (p=7.6e-05)
   against the regression's 0.038 (p=0.41), but still roughly a third of the
   other two regions. This does **not** rehabilitate the caudate cell. Matching
   on GC and coverage cannot fix a variable that is perfectly confounded with
   region any more than adjusting for it can; it only removes the extrapolation.
   The caudate estimate stays labelled technically uncertain and stays
   non-gating. What it does establish is that the caudate null was partly an
   artifact of the adjustment's functional form, so "caudate shows nothing" is
   not a claim these data support either.

## QC scripts

`_h/06` and `_h/07` are post-hoc analyses OF sealed runs, not stages of one.
They write to `_m/qc/{run_id}/`, which is gitignored and regenerable, and
neither can produce a number that enters `association-results.tsv` or the gates.
`_h/08` is different: it runs as part of step 2 and writes into the run's
`results/`, but it is still non-gating by construction -- its output files are
never read by `03_apply_gates.R`.

| script | question |
|---|---|
| `_h/06_qc_adjustment_ladder.R` | Which BLOCK of covariates carries the v1 -> v2 attenuation? Cumulative nested ladder, plus the prespecified set fitted absolutely. |
| `_h/07_qc_covariate_attribution.R` | Which single covariate, and does it act the same way in every region? Four passes: `attribution` (add-one / leave-one-out), `correlation` (predictor and outcome vs every covariate), `coverage` (minimum-depth sweep), `set_overlap` (shared vs region-unique VMRs). `--outcome` defaults to `line_l1_frac`. |

| `_h/08_matched_measurability.R` | Does the association survive replacing regression adjustment with matching on measurability? Non-gating; part of step 2, writes `descriptive-matched-measurability.tsv` and `descriptive-matched-balance.tsv` into the run. |

Read `06` and `07` as decompositions, never as a menu to select an adjustment
set from.

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
