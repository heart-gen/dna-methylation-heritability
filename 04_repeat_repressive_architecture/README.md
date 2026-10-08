# 04_repeat_repressive_architecture — primary biological analysis

Tests whether a higher relative local SNP contribution score (`local_snp_contribution_score_z`, Module 02) is associated with repeat-rich and repressive genomic compartments. This is the module the manuscript's central claim rests on.

**Status: accepted (AA, 2026-09-25), claims as in the Accepted runs table.** The
accepted runs are `rra-AA-{caudate,dlpfc,hippocampus}-20260925-a`. Both corrections
recorded under **Corrections landed 2026-09-23** are in effect in them, they carry
the BrainScope ATAC contrast added the same day, and every earlier acceptance is
superseded.

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
`00_new_run.R` with the message. The smoke runs above were created
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

## Corrections landed 2026-09-23 (a rerun is required for either to take effect)

Both were found in the PI-summary pass, both change sealed outputs, and neither
is reflected in the accepted `rra-AA-*-20260906` runs.

### 1. The composition adjustment is RNA MuSiC everywhere; scMD is caudate-only

`_h/01_build_features.R` loaded `dnam-scmd-proportions-{region}.tsv`
unconditionally and no script in the module read RNA MuSiC at all, so the single
`cell_composition_r2` covariate -- and both gating composition arms built on it
-- were scMD-derived in all three regions. AGENTS.md 7.4 is asymmetric: the RNA
MuSiC adjustment is required unconditionally, the DNAm scMD adjustment only
"when the integration gate passes". That gate passes in caudate
(neuronal rho 0.720, FDR 2e-45) and fails in DLPFC (-0.035, p 0.66) and
hippocampus (-0.103, p 0.27), so two of three regions were adjusted for a
deconvolution that does not track the composition it claims to measure.
`config/repeat_annotations.yml:sensitivities` already declared both
`cell_composition_rna_music: true` and
`cell_composition_dnam_scmd: true  # caudate only, when the integration gate
passes`; the locked configuration asked for this before the code did it, and no
config change was needed to implement it.

What changed:

| column / arm | before | after |
|---|---|---|
| `cell_composition_r2` | DNAm scMD, all regions | **RNA MuSiC, all regions** |
| `cell_composition_r2_music` | -- | RNA MuSiC (same values, named) |
| `cell_composition_r2_scmd` | -- | DNAm scMD where the gate passes, `NA` elsewhere |
| `cell_composition_r2_source`, `scmd_integration_gate` | -- | provenance, on every row |
| `low_cell_composition`, `adjust_cell_composition` | scMD-derived, gating | **MuSiC-derived**, gating in all three regions |
| `adjust_cell_composition_scmd` | -- | new arm, **fitted in caudate only**, gating where fitted |
| `results/cell-composition-arms.tsv` | -- | which modality was used, from which file, and why not |

`cell_composition_r2` keeps its name because PI-locked configuration downstream
refers to it (`config/aging.yml:axis.sensitivities`,
`config/gwas_negative_controls.yml:continuous_z`); its modality is now recorded
in the table rather than implied by the script. **Those downstream arms change
value**, so 09b and Module 09 stage 18 must be rerun after this module.

Where the gate fails, the scMD arm is present in `association-results.tsv` with
`arm_fitted = FALSE` and the reason `scmd_integration_gate_fails_in_region`,
never silently absent. `03_apply_gates.R` excludes unfitted arms from the
survival conjunction -- counting them would fail every outcome in two regions --
and records which arm was missing where in the claims table's
`sensitivity_arms_not_fitted`. `08_region_donor_generalization` drops unfitted
rows at harvest, so a sensitivity that was never run is not reported as one that
did not replicate. The gate token's denominator is unaffected: it counts
outcomes (3), not arms.

Shared code: `00_shared/cell_composition.R` now holds `scmd_gate_passes()`
(moved unchanged from `09b_aging_application/_h/age_functions.R`, which defined
it first), plus `cell_composition_sources()` and `vmr_composition_r2()`.

#### Closed 2026-09-27: the scMD arm gates, and on the accepted run it costs nothing

`cell_composition_dnam_scmd` carries no `gating:` key, and in this config block
every sensitivity gates unless it says otherwise -- `exclude_snp_proximal_cpgs`
has `descriptive_only: true` and `matched_measurability` has `gating: false`. The
implementation follows that reading, so **caudate's survival conjunction carries
one arm that DLPFC's and hippocampus's do not.** The asymmetry is real and is
named per outcome in the claims table's `sensitivity_arms_not_fitted`.

The question was whether to add `gating: false` so a caudate-only arm can never
break a shared claim. It was left as written, because on
`rra-AA-*-20260925-a` the arm changes no verdict, and that is checkable rather
than asserted:

| caudate outcome | scMD-adjusted arm | verdict, and what actually decides it |
|---|---|---|
| `quiescent_frac` | est 0.340, p = 6.5e-28 | survives; **3/3 supported** and the scMD arm is not the weakest link |
| `h3k9me3_frac` | est 0.093, p = 0.116 | caudate already fails at `high_mappability` (p = 0.545) and `adjust_cell_composition` (p = 0.141), so **2/3** holds with or without the scMD arm |
| `line_l1_frac` | -- | caudate is set aside as technically confounded (§8.1) and is not counted at all |

So the gating decision is presently **moot**: adding `gating: false` would change
no claim, no token and no figure, and would change the config checksum. It stays
as locked. This becomes live again the moment a caudate outcome's only failing arm
is the scMD one, so the check above should be rerun whenever the claims table
moves.

One interaction to carry forward, because it is not visible from this module: the
paragraph above says Module 08 drops unfitted rows at harvest so that "a
sensitivity that was never run is not reported as one that did not replicate."
That is true row by row and **not** true of the conjunction built on top of them.
Dropping DLPFC and hippocampus leaves the scMD arm with `n_regions = 1`, and
`08/_h/01_cross_region_replication.R` requires all three regions for
`complete_across_regions`, so the arm can never satisfy strict replication and
every repeat test inherits `replicated_strict = FALSE`. See
`08_region_donor_generalization/README.md`, "The accepted run now rests on
superseded upstreams" -- the fix belongs there, not here.

### 2. The run decision token counted prose, not gates

`_h/05_finalize_run.R` derived its support count as
`sum(!startsWith(claims$permitted_claim, "not supported"))`. H3K9me3's claim
begins "below the gate (2/3 regions)", so it was counted as supported and all
three cells sealed as `GATES_APPLIED_3_OF_3_OUTCOMES_SUPPORTED` when **two**
outcomes had cleared. The claims table, this README and
`MIGRATION_MANIFEST.tsv` recorded the correct 2-of-3 throughout; only the token
was wrong, and it was wrong because it parsed a sentence.

`03_apply_gates.R` now emits the structured verdict -- `gate_supported`
(`regions_surviving >= regions_required`) and `gate_status` (`supported` /
`below_gate` / `not_supported`) -- validates the gate counts, cross-checks them
against the claim text, and `05_finalize_run.R` counts that column and refuses
to seal a claims table that lacks it. The sealed 2026-09-06 manifests keep their
overstated token; read `interpretation-claims.tsv`, not the token, until the
rerun.

## Shared vs region-unique VMR split (non-gating, 2026-10-07)

Run `rra-AA-crossregion-20261007` was produced by `_h/09_shared_unique_split.R`
at commit `a233c4543`. It reads the three accepted `-20260925-a` cells, and its
rebuilt primary model reproduces all 45 sealed primary fits to 1e-6. A VMR is
**shared** when it overlaps (≥ 1 bp) a VMR in each of the other two regions.
Every other VMR is **region-unique**.

| | caudate | DLPFC | hippocampus |
|---|---:|---:|---:|
| shared VMRs | 3,370 | 3,409 | 3,417 |
| region-unique VMRs | 7,881 | 5,842 | 5,749 |

Values are the primary-model coefficient on `local_snp_contribution_score_z`,
with the chromosome-jackknife p in parentheses.

| outcome | subset | caudate | DLPFC | hippocampus |
|---|---|---:|---:|---:|
| quiescent | shared | +0.303 (5.3e-05) | +0.246 (1.0e-04) | +0.292 (4.2e-07) |
| quiescent | unique | +0.541 (2e-48) | +0.501 (1e-54) | +0.570 (6e-55) |
| H3K9me3 | shared | −0.140 (0.32) | −0.035 (0.77) | −0.058 (0.54) |
| H3K9me3 | unique | +0.253 (1.9e-04) | +0.543 (2.8e-12) | +0.436 (1.2e-09) |
| LINE/L1 | shared | +0.106 (0.20)† | +0.060 (0.49) | +0.093 (0.19) |
| LINE/L1 | unique | −0.026 (0.71)† | +0.430 (3.4e-15) | +0.455 (3.3e-16) |
| accessible | shared | −0.289 (4e-11) | −0.189 (1.2e-04) | −0.234 (7.0e-06) |
| accessible | unique | −0.559 (2e-72) | −0.606 (4e-63) | −0.631 (2e-56) |
| ATAC union | shared | −0.264 (1.3e-12) | −0.304 (1.3e-10) | −0.255 (8.3e-09) |
| H3K27ac | shared | −0.205 (3.2e-06) | −0.246 (7.8e-11) | −0.250 (3.8e-07) |

† Caudate LINE/L1 stays set aside from the claim, as in the primary.

The joint-jackknife unique-minus-shared difference for H3K9me3 is +0.39 /
+0.58 / +0.49 (p 0.016 / 3.4e-05 / 1.0e-07). For LINE/L1 in DLPFC and
hippocampus it is +0.37 / +0.36 (p 7.7e-06 / 2.4e-05).

**Reading.**
- **Present in shared VMRs in all three regions:** the quiescent enrichment and
  the depletion from accessible chromatin (DNase, ATAC union, H3K27ac). Their
  strength there is about half of what it is in region-unique VMRs. These are
  properties of the VMR set as a whole.
- **Absent from shared VMRs in every region:** the H3K9me3 association and the
  DLPFC and hippocampus LINE/L1 associations. Region-unique VMRs carry all
  three, and the gap survives the joint jackknife.

This run cannot separate two readings:
1. Genetically controlled LINE/L1 and H3K9me3 methylation is genuinely
   region-restricted.
2. Region-unique VMRs are enriched for loci whose calling, coverage or score is
   sensitive to local sequence.

Neither reading changes a gate. The Figure 3 sentence should not describe the
LINE/L1 or H3K9me3 association as a property of VMRs common to all three
regions.

## LINE/L1 subfamily resolution (non-gating, 2026-10-07)

**Why.** §7.4 lists, as a high-value extension, separating LINE/L1 by subfamily,
comparing young L1HS/L1PA with older L1M, and distinguishing full-length from
truncated elements. `repeat_annotations.yml:repeatmasker.l1_subfamily_definitions`
has been `null` since 2026-08-23. The PI asked for it on 2026-10-07.

**Definitions** are in `config/l1_subfamilies.yml` (`pi_locked`). They were
fixed before any subfamily association was fitted; only element counts were
looked at. The file is separate so that `repeat_annotations.yml` keeps the hash
the accepted cells recorded.
- **Age**, by RepeatMasker subfamily name, first match wins:
  - young `^(L1HS|L1PA[0-9])`;
  - old `^L1M`;
  - intermediate, every other `^L1` (L1P1-5, L1PB*, L1PREC2). It is fitted
    but kept out of the young-vs-old contrast.
- **Full-length:** the subfamily consensus is ≥ 5.5 kb, and the element's
  alignment covers ≥ 90% of it.
- **Retains 5′ end:** the alignment starts within the first 100 bp of the
  consensus. Descriptive only.
- **Contrasts**, each a difference of slopes on the same VMRs with a joint
  chromosome jackknife: young − old, full-length − fragment, and young
  full-length − young fragment. The last separates completeness from age,
  because 8,249 of the 9,247 full-length elements are young.

**Asset.** `inputs/supportfiles/_h/02_build_l1_subfamily_asset.py` reads UCSC
hg38 `rmsk.txt.gz` (sha256 pinned in the config) and writes
`repeat-masker.LINE_L1.subfamily.hg38.tsv.gz`: 989,683 elements in 127
subfamilies. It first proves the UCSC table is the one the project assets came
from, by sorted MD5 of the six-column projection against
`repeat-masker-hg38.gz` and of the LINE/L1 subset against
`repeat-masker.LINE_L1.hg38.bed.gz`. Both match
(`_m/l1-subfamily-asset-report.tsv`). Consensus coordinates are
strand-dependent in the UCSC table, and the builder reads them that way.

| | elements | full-length | fragment |
|---|---:|---:|---:|
| young | 130,658 | 8,249 | 122,409 |
| intermediate | 57,051 | 558 | 56,493 |
| old | 801,974 | 440 | 801,534 |

**Stage.** `_h/10_l1_subfamily.R --cohort AA` reads the three accepted cells
and fits each class fraction with the primary model in the four analysis
sets. Two guards must pass before anything is written:
- **Partition guard:** every element falls in exactly one age class, and for
  every VMR the union of the classes, recomputed from the subfamily asset, must
  reproduce the sealed `line_l1_frac` to 1e-12. The class fractions themselves
  can sum to slightly more than `line_l1_frac` where elements of two classes
  overlap, as with a young L1 inserted into an old one (up to 0.08 of a VMR's
  span; `l1-subfamily-partition-check.tsv`). The config comment that says they
  sum exactly is imprecise; the guard is on the union.
- **Reproduction guard:** the rebuilt model must reproduce the sealed
  `line_l1_frac` primary fit to 1e-6.

An outcome is fitted in an analysis set only if at least 100 VMRs overlap the
class there. Otherwise the row is written as `insufficient_overlap` with no
estimate. The run is `rra-{cohort}-crossregion-{date}`, the same token as the
shared/unique split; the manifest's `sensitivity =
line_l1_subfamily_resolution` tells them apart.

**Reading rule**, fixed before fitting. A contrast is read as `consistent` only
if all three conditions hold:
- in the primary set, its joint-jackknife p < 0.05 in both DLPFC and
  hippocampus;
- the sign is the same in both regions;
- the high-mappability estimate keeps that sign in both.

Caudate is fitted and shown but set aside, as caudate `line_l1_frac` is (§8.1).
A contrast whose high-mappability arm falls under the overlap floor reads
`not_evaluable_high_mappability`, never `consistent`.

**Limits.**
- Overlap does not show activity, expression or retrotransposition.
- Full-length is not retrotransposition competence, since there is no ORF
  integrity check; `retrotransposition_competent_list` stays `null`.
- The UCSC table does not link the two pieces of an element split by an
  insertion, so the full-length class is conservative.
- Young L1 is the least mappable sequence in the genome, which is why the
  high-mappability arm is in the reading rule.

**Result: `rra-AA-crossregion-20261007-a`** (sealed 2026-10-07 at `b924291e3`,
clean tree, host run). The union check passed to 5.6e-16 and the model
reproduction to ≤ 4.4e-16 in all three regions. All three contrasts read
`not_consistent`.

Values are the coefficient on `local_snp_contribution_score_z`, with the
chromosome-jackknife p and the number of VMRs overlapping the class.

| outcome | DLPFC | hippocampus | caudate (set aside) |
|---|---|---|---|
| old L1M | +0.372 (7e-16; 954) | +0.403 (3e-11; 997) | +0.314 (1e-04; 1,224) |
| young L1HS/L1PA | +0.285 (0.008; 360) | +0.296 (0.006; 357) | −0.251 (0.003; 497) |
| intermediate | +0.207 (0.3; 122) | +0.122 (0.5; 119) | +0.284 (0.2; 157) |
| fragment | +0.334 (5e-13; 1,239) | +0.363 (3e-17; 1,279) | +0.121 (0.04; 1,519) |
| full-length | +0.316 (0.1; 114) | +0.236 (0.2; 113) | −0.084 (0.4; 183) |
| retains 5′ end | +0.504 (0.01; 155) | +0.523 (1e-04; 155) | +0.111 (0.4; 228) |
| old L1M, high mappability | +0.322 (3e-07; 482) | +0.357 (1e-06; 504) | +0.197 (0.04; 640) |
| young, high mappability | not fitted (56) | not fitted (56) | not fitted (98) |

| contrast (joint jackknife) | DLPFC | hippocampus | caudate (set aside) |
|---|---|---|---|
| young − old | −0.086 (p 0.40) | −0.107 (p 0.36) | −0.564 (p 1.1e-05) |
| full-length − fragment | −0.019 (p 0.93) | −0.126 (p 0.49) | −0.205 (p 0.12) |
| young full-length − young fragment | +0.095 (p 0.75) | −0.064 (p 0.79) | −0.064 (p 0.71) |

**Reading.**
- **In DLPFC and hippocampus the LINE/L1 association is not specific to young
  or to full-length elements.** Old L1M and young L1HS/L1PA both rise with the
  score, at slopes that do not differ. Old L1M overlaps about 2.7 times as many
  VMRs and carries most of the precision. It is also the only class that can
  be tested under high mappability, and it survives there.
- **Young L1 cannot be tested for mappability robustness.** Only 56 VMRs per
  region overlap young L1 after the high-mappability restriction, under the
  100 floor. Its primary estimate therefore stays unverified against the
  mapping artefact the restriction exists to catch. Full-length and
  5′-retaining classes overlap 8-53 VMRs there in DLPFC and hippocampus and are
  untestable too.
- **Full-length versus fragment is underpowered, not null.** About 114 VMRs
  overlap a full-length element in each region, and the full-length slope's
  jackknife SE is 0.2. The 5′-retaining class is nominally significant in
  both regions but is not a prespecified contrast.
- **Caudate is set aside, but its split is informative about the collapse.**
  Caudate `line_l1_frac` is null (+0.016) because a positive old-L1M slope
  (+0.314) and a negative young-L1 slope (−0.251) cancel. DLPFC and hippocampus
  show no such sign split. Caudate is sequencing batch 3 (§8.1), and young L1 is
  the least mappable class, so this is consistent with the batch acting through
  mappability-sensitive sequence. It cannot be read as regional biology.
- Nothing here licenses a statement about activity, expression or
  retrotransposition. The main-text LINE/L1 sentence is unchanged. At most it
  can add that the association is carried mainly by old L1M sequence and does
  not distinguish young from old elements.

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
| `_h/09_shared_unique_split.R` | Is each primary association a property of the region's whole VMR set, or carried by the VMRs only that region called? `07`'s `set_overlap` split, made citable: it reads the three accepted cells and refuses to fit unless it reproduces every sealed primary estimate to 1e-6. It then seals its own cross-region run, `rra-{cohort}-crossregion-{date}`. Each subset gets HC3 and chromosome-jackknife SEs; the unique-minus-shared difference uses a joint jackknife. Non-gating: it changes no claim, q or decision token. |
| `_h/10_l1_subfamily.R` | Does the LINE/L1 association differ between young L1HS/L1PA and old L1M elements, or between full-length elements and fragments? Reads the three accepted cells, refuses to fit unless the age classes partition `line_l1_frac` (1e-12) and the rebuilt model reproduces its sealed fit (1e-6), then seals `rra-{cohort}-crossregion-{date}`. Non-gating and outside the BH family; see the section above. |

Read `06` and `07` as decompositions, never as a menu to select an adjustment
set from.

`06` and `07` select production cells by the pattern
`rra-AA-(caudate|dlpfc|hippocampus)-YYYYMMDD`. That pattern excludes the
cross-region runs (the split and the L1 subfamily stage), and it also excludes the accepted `-20260925-a` cells.
Pass `--run-id` to `06` to target an accepted cell.

## Contract

This module follows: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.

## Accepted runs

The three cells were gated together; `interpretation-claims.tsv` lives on the
caudate run. Both 2026-09-23 corrections are in effect here, and each is visible
in the sealed output:

- the composition adjustment is RNA MuSiC in all three regions, with the DNAm scMD
  arm fitted in caudate only -- `sensitivity_arms_not_fitted` reads
  `adjust_cell_composition_scmd not fitted in dlpfc,hippocampus
  (scmd_integration_gate_fails_in_region)` rather than leaving the arm silently
  absent;
- the decision token is derived from the structured `gate_supported` column rather
  than parsed from claim prose, so it reads
  `GATES_APPLIED_2_OF_3_OUTCOMES_SUPPORTED` where the 2026-09-08 manifests
  overstated 3 of 3. The claims table, this README and `MIGRATION_MANIFEST.tsv`
  recorded the correct 2 of 3 throughout; only the token was wrong.

Claims licensed, unchanged in substance from 2026-09-08. Written as a list rather
than a table on purpose: `00_shared/gates.R::read_accepted_runs()` takes the
**first** markdown table under this heading as the acceptance record and keeps only
rows matching that table's column count, so a second table here silently replaces
the record with itself.

- `quiescent_frac` -- gate 3 of 3 -- **supported**, shared across all three regions.
- `line_l1_frac` -- gate 2 of 2 eligible -- **supported** in DLPFC and
  hippocampus; caudate set aside as technically confounded (estimate 0.0161,
  p = 0.729), reported and not counted.
- `h3k9me3_frac` -- gate 3 of 3 -- **below the gate** (2 of 3); describable only as
  suggestive in DLPFC and hippocampus.

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| rra-AA-caudate-20260925-a | AA | caudate | vmrset-AA-caudate-937a41979978 | 2026-09-25 | Kynon J.M. Benjamin | GATES_APPLIED_2_OF_3_OUTCOMES_SUPPORTED | Gate host; `interpretation-claims.tsv` and `association-results-all-regions.tsv` live here. Quiescent 3/3; LINE/L1 caudate excluded from the claim (0.0161, p=0.729); H3K9me3 below the shared gate. Caudate remains GC-entangled. 16 ATAC outcomes added outside the family |
| rra-AA-dlpfc-20260925-a | AA | dlpfc | vmrset-AA-dlpfc-856067dfe289 | 2026-09-25 | Kynon J.M. Benjamin | GATES_APPLIED_2_OF_3_OUTCOMES_SUPPORTED | Survives quiescent and LINE/L1; H3K9me3 suggestive only. scMD arm not fitted (integration gate fails) |
| rra-AA-hippocampus-20260925-a | AA | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-09-25 | Kynon J.M. Benjamin | GATES_APPLIED_2_OF_3_OUTCOMES_SUPPORTED | Survives quiescent and LINE/L1; H3K9me3 suggestive only. scMD arm not fitted (integration gate fails) |

Provenance: `vmrcat-AA-{region}-20260816` -> `lgv-AA-{region}-rescore-20260913` ->
`lsp-AA-{region}-20260925-a` -> this run, sealed 2026-09-25T23:09 at commit
`d61f83b4c`, `git_dirty = false`, `smoke_run = FALSE`.

### Accepted non-gating sensitivity runs

These runs read the three accepted cells above and change no claim, q-value or
decision token. They are recorded here, under their own heading, because
`read_accepted_runs()` allows one accepted run per cohort x region and both are
`AA x crossregion`. Nothing downstream gates on them, and Module 11's
analysis-to-claim matrix does not list them. Each run's manifest field
`sensitivity` says which analysis it is.

| run_id | sensitivity | upstream | sealed | accepted_on | accepted_by | notes |
|---|---|---|---|---|---|---|
| rra-AA-crossregion-20261007 | shared_vs_region_unique_vmrs | rra-AA-{caudate,dlpfc,hippocampus}-20260925-a | 2026-10-07T14:58 at a233c4543, git_dirty false | 2026-10-08 | Kynon J.M. Benjamin | Reproduces all 45 sealed primary fits to 1e-6. Quiescent enrichment and accessible-chromatin depletion are present in shared VMRs in all three regions, at about half the strength they have in region-unique VMRs. The H3K9me3 and the DLPFC/hippocampus LINE/L1 associations are absent from shared VMRs and carried by region-unique ones; the unique minus shared gap survives the joint jackknife (H3K9me3 p 0.016 / 3.4e-05 / 1.0e-07; LINE/L1 p 7.7e-06 / 2.4e-05). The run cannot separate region-restricted biology from region-unique VMRs being more sequence-sensitive. Figure 3 must not call LINE/L1 or H3K9me3 a property of VMRs common to all regions. |

### The ATAC contrast, added 2026-09-25 (T8 and T11)

These runs are the first to carry the BrainScope ATAC CRE tracks: one published
union (`atac_union_frac`, role `complementary_contrast_independent_assay`) and seven
per-cell-type tracks (`atac_{astro,endo,exc,inh,micro,opc,oligo}_frac`, role
`celltype_breakdown_secondary`). All sixteen `_frac`/`_any` outcomes sit in
`outside_family`, gated on a one-sided raw `p` against a declared negative
direction, and none of them touches a reported q. That is verified rather than
assumed: the nine BH-family rows are **bit-identical** to the superseded
`-20260925` run in estimate, `p` and `q`, while the cross-region results file grew
from 504 to 1,080 rows. Adding a control cannot revise a claim, and here it
demonstrably did not.

**What the union buys.** `accessible_frac` is Roadmap ChromHMM; `atac_union_frac`
is ATAC-seq from a different consortium, different donors and a different assay
chemistry. They agree closely:

| region | `accessible_frac` | `atac_union_frac` |
|---|---|---|
| caudate | -0.437 | -0.387 |
| DLPFC | -0.431 | -0.496 |
| hippocampus | -0.449 | -0.433 |

So the depletion of local genetic control in accessible chromatin is not a property
of one annotation pipeline. This is the one genuinely new piece of evidence in the
rerun, and it strengthens a **control**, not a claim: the accessible-chromatin
contrast qualifies the repressive-compartment result by showing the gradient runs
the other way in active sequence, which is what a complementary contrast is for.

**What the seven cell types do not buy.** Every one is negative and passes its
one-sided gate in all three regions, with magnitudes from -0.19 (microglia,
hippocampus) to -0.44 (OPC and astrocyte). **No cell type separates from the
others**, and the breakdown must not be read as identifying one:

- §2.3 forbids inferring a cell type of origin from bulk tissue, and nothing about
  this design escapes that. The VMRs are bulk WGBS; only the annotation is
  cell-resolved.
- The seven tracks are not independent. Pairwise Jaccard is 0.16-0.33, so a
  consistent sign across all seven is close to the expected outcome for any
  genome-wide accessibility gradient, not seven concurring tests.
- The seven do not reconstruct the union. Only 93.6-97.9% of each cell type's
  intervals overlap the published union track, so the union is a separate
  observation rather than a summary of the breakdown -- which is why both are
  registered.

Read the seven as a check that the union result is not driven by a single cell
type's peaks, and stop there.

### Resolved 2026-09-25: the secondary predictor is now on the rescored Module 03

`upstream_local_snp_prediction_run_id` reads `lsp-AA-{region}-20260925-a` in all
three cells. The earlier acceptances pointed at `lsp-AA-{region}-20260825`, which
predated the Module 03 rescore reruns, so the primary predictor was on the rescore
and the secondary `r2_pred_oof_z` arm was not.

That mismatch is closed, and the check that was outstanding has now been run. It
never threatened a licensed claim -- the gate already refuses to let that arm carry
weight, and the LINE/L1 claim reads, in the sealed claims table: *"rests on the
primary predictor ALONE -- the `r2_pred_oof_z` association is descriptive and
near-circular, and may not be cited as corroboration."* The secondary arm lives in
its own file (`secondary-predictor-descriptive.tsv`) and never enters the survival
conjunction. What it did mean was that no number in that file could be quoted
without a refresh first. Numbers from `secondary-predictor-descriptive.tsv` in
`-20260925-a` are quotable; numbers from any earlier run are not.

### Superseded

`rra-AA-{caudate,dlpfc,hippocampus}-20260925` (accepted 2026-09-25, superseded the
same day). Sound runs, superseded only because they predate two things: the
BrainScope ATAC registration, and the Module 03 rescore pointer. Their BH family is
bit-identical to the accepted runs, so every claim they licensed is unchanged --
they are superseded for completeness of the outcome set and provenance, not for
error. Do not cite them; nothing in them is wrong.

`rra-AA-{caudate,dlpfc,hippocampus}-20260906` (accepted 2026-09-08). Superseded
2026-09-25 on three counts, none of which changed a licensed claim: the scMD/MuSiC
composition defect (two of three regions were adjusted for a deconvolution that
does not track the composition it claims to measure), the prose-parsed decision
token, and the Module 02 pointer (`lgv-AA-{region}-20260823`, replaced by the
rescore). Their `interpretation-claims.tsv` remains readable for audit; their
manifests' `GATES_APPLIED_3_OF_3_OUTCOMES_SUPPORTED` token must not be cited.

`rra-AA-{caudate,dlpfc,hippocampus}-20260902` -- the first production run, never
accepted. Superseded by `-20260906` after the 2026-09-02 adjustment-set amendment
and the `exclude_snp_proximal` retirement. Still referenced above by the `_h/06`,
`_h/07` and `_h/08` QC analyses, which were run against its sealed feature tables.

`rra-smoke-AA-{region}-20260823` and `rra-smoke2-AA-{region}-20260925` -- smoke
runs, `smoke_run = TRUE`, sealed with
`GATES_APPLIED_..._SMOKE_ONLY_NOT_ACCEPTABLE`. Never citable.
