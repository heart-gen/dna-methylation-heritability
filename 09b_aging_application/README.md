# 09b_aging_application — age-associated methylation along the local-genetic-control axis

**Status: implemented and smoke-tested 2026-09-19. `config/aging.yml` locked
by the PI 2026-09-19; production runs `age-AA-{region}-20260919` accepted the
same day (see Accepted runs). Cross-region token: NOT_SUPPORTED.**

## Question

Do VMRs with **weaker** local SNP control show **larger** age-associated
methylation differences?

This is an orthogonal, non-disease test of the pattern Module 09 found for
schizophrenia, where SCZ-linked VMRs sit lower on the local-genetic-control axis
in all three regions. If lower genetic control marks methylation that responds
more to aging, the axis carries biology beyond disease. The module is kept
lean on purpose: one axis test per region plus one cross-region stage. It is
not a second disease module.

The `09b` prefix follows the `01b` precedent. The module depends on 01, 02, 04
and 07 and **not** on 09; it sits beside 09 as the second application of the
axis. It is not a sub-stage of 09. Module 07 is read only by the descriptive
annotation table.

## Design

### Per-VMR age model (stage 01)

`meth ~ age + sex + diagnosis` is fitted for every Module 02-eligible VMR, with
age in decades. The covariates come from
`00_shared/locus_io.R::load_locus_phenotype()`, the reader Module 02 used, so
the age effect conditions on exactly what the score conditions on.

| spec | model | role in the region reading |
|---|---|---|
| `primary` | age + sex + diagnosis | primary |
| `controls_only` | age + sex, Control donors | gating (reduced n) |
| `cell_music` | primary + RNA MuSiC cellPC1–3 | gating |
| `cell_scmd` | primary + DNAm scMD cellPC1–3 | gating, **caudate only** (scMD passes its integration gate only there) |
| `age_ge_25` | primary, donors ≥ 25 years | non-gating |

Cases are about 7 years older than controls. The AA medians are 53.7 vs 48.6
(caudate), 54.3 vs 45.8 (hippocampus) and 53.2 vs 45.8 (DLPFC), over a range of
17–84. Diagnosis is therefore a covariate, and every bootstrap draw resamples
within diagnosis.

Cell PCs come from `00_shared/cell_composition.R::cell_composition_pcs()`. It
was lifted out of Module 04's feature builder so both modules construct them
identically; it reproduces Module 04's PCs bit for bit.

### The outcome, and why it is not |β| (stage 02)

The primary axis model is

    (β̂² − SE²) / mean  ~  local_snp_contribution_score_z + vmr_length + cpg_count
                          + cpg_density + gc_content + mappability + mean_methylation

The per-VMR age model has no SNP term. A VMR with strong local SNP control
therefore carries that genetic variance in its residual and gets a larger SE.
Any age-effect measure whose expectation depends on the SE is coupled to the
score mechanically:

- **|β̂| and its rank.** E|β̂| ≈ √(β² + SE²), so the noise floor rises with the
  score.
- **|t|, −log10 p and partial R².** These are pivotal when there is no age
  effect. Once real effects exist, a fixed β gives a smaller t where the SE is
  larger.

**β̂² − SE² has expectation β² whatever the SE**, so the coupling is removed by
construction rather than modelled. Dividing by the region mean makes the
coefficient read as "proportional change in mean squared age effect per SD of
score". That quantity is unit-free, which is what makes the DLPFC–hippocampus
comparison meaningful.

`methylation_variance` is kept out of the primary covariates. Total variance
contains the age variance, so adjusting for it would partly adjust away the
outcome. It returns as a non-gating arm.

**Inference.** SE² = donor-bootstrap variance (within diagnosis) + weighted
delete-one-chromosome block-jackknife variance (Busing et al. 1999, as in
Modules 06 and 08).

- The donor half accounts for shared donors correlating every VMR's estimate.
- The VMR half accounts for the VMRs being a sample of correlated loci.
- Either half alone under-covers.
- The bootstrap is used for its variance only. Its distribution of β̂² − SE² is
  shifted by about SE², so a percentile interval would reintroduce the coupling.

### Rejected on the same day: rank(|β|) + age permutation

The first design, approved in planning, used the rank of |β̂| with a
donor-level permutation of age. The observed coefficient was read against the
permutation null.

Simulation showed it is invalid:
- The permutation null sets **every** age effect to zero, which is where the
  noise-floor coupling is strongest.
- With real age effects **unrelated** to the score, the observed coefficient
  sat below the null mean every time.
- It rejected **100%** of such nulls, in the hypothesized direction.

A permutation of age tests "no age effect anywhere", which is not this
module's null. The simulations used real DLPFC donors and 1,500 VMRs with
genetic SD rising with score. The simulation scripts are not tracked; the
retained checks are in `tests/`:

| condition | rank(\|β\|) + permutation | debiased, donor bootstrap only | debiased, combined SE |
|---|---|---|---|
| effects unrelated to score, τ = 0.01 | 100% reject | 15% | 0% |
| effects unrelated to score, τ = 0.03 | 100% reject | 57% | 3% (5% in `tests/`) |
| effects shrink with score (grad 0.25–0.5) | — | 100% | 100% |

### Arms and secondaries (stage 02)

All of these are fitted on the primary spec's age effects.

| arm / outcome | role |
|---|---|
| + `cell_composition_r2` | gating; **caudate only when the column is scMD-derived**, every region when it is MuSiC-derived. The modality is read from Module 04's `cell_composition_r2_source`, not from the name (see below) |
| + `methylation_variance` | non-gating |
| mappability ≥ 0.9 | non-gating |
| + H3K27me3, bivalent, H3K9me3 and quiescent fractions | non-gating decomposition: does the gradient ride on PRC2/bivalent age-hypermethylation or on Module 04's architecture? Attenuation links the two; it is not a failure. |
| signed β̂ (gain vs loss) | secondary, same inference; β̂ is unbiased, so no debiasing is needed |

### Annotation associations (stage 02, descriptive)

Which kinds of VMR carry the age-associated differences? This is the other half
of the question the module asks. The Module 04 architecture and the Module 07
coupling are what distinguish low-control from high-control VMRs, so how age
effects distribute over them is biology in its own right. It is not treated as
a nuisance to adjust away.

For each annotation, one at a time, the same two outcomes are regressed on
that annotation instead of the score:
- magnitude: relative β̂² − SE²;
- direction: signed β̂.

The models use the primary spec's age effects, the same donor-bootstrap draws
and the same combined SE.

| source | annotations |
|---|---|
| Module 04 chromatin (any overlap) | H3K27me3, bivalent, H3K9me3, quiescent, accessible, H3K27ac |
| Module 04 repeats | LINE/L1 |
| Module 04 genomic context | promoter, exonic, intronic, intergenic (each vs the rest) |
| Module 04 continuous | `cell_composition_r2` (z) |
| Module 07 | coupled to nearest-gene expression, PSI, ABC-linked expression (adjusted for features tested; ≥ 50 coupled VMRs, otherwise recorded as not tested) |

How to read the table:
- **Two adjustments.** Every annotation is fitted with the technical axis
  covariates, and again with the score added. The second fit asks whether the
  annotation's age association is carried by local genetic control or is
  separate from it.
- **Multiple testing.** BH q is taken within region × outcome × adjustment.
- **No mutual adjustment.** Annotations are not adjusted for one another, and
  promoter, H3K27ac and accessible overlap.
- **Unadjusted contrasts.** The unadjusted means and the fraction of VMRs that
  gain methylation are reported beside each estimate.
- **Cross-region.** Stage 05 counts the regions with q < 0.05 outside caudate
  and records whether their signs agree. Caudate is shown beside them, not
  counted.

**Provenance.** This table was prespecified on 2026-09-19, *after* an
exploratory pass over the smoke runs' per-VMR age effects against these
annotations. That pass did not look at the score–age relation.
- To limit selection, the list is every biological annotation class Modules 04
  and 07 emit, not a subset chosen on the results.
- Technical flags (segdup, problematic, SNP-proximal) are excluded.
- The exploratory pass found:
  - H3K27me3, bivalent, promoter and H3K9me3 VMRs **gain** methylation with
    age in all three regions;
  - accessible and H3K27ac VMRs **lose** methylation;
  - intergenic VMRs have smaller age effects;
  - Module 07 coupling showed no association.

  Treat these as the reason the table exists, not as its result.

### Region reading (stage 03)

The primary must be significant (p < 0.05) with a **negative** coefficient.
Then a strict conjunction applies:
- `cell_music`, `cell_scmd` and `cell_composition_r2` have the same n, so each
  must keep the sign **and** stay significant.
- `controls_only` keeps about 58% of donors, so it must keep the sign and at
  least 50% of the primary coefficient.

A gating member that cannot be fitted in a region is recorded as absent
evidence, not contrary evidence, with its reason in
`gating-sensitivities.tsv:reason`.

**Corrected 2026-09-23.** A gating arm whose covariate is **scMD-derived** is
declined where scMD fails its integration gate, exactly as the `cell_scmd` spec
is, and recorded in `axis-arms-skipped.tsv`. Previously `cell_composition_r2` —
then the per-VMR R² of methylation on the *scMD* proportion PCs, the same
quantity `cell_scmd` adjusts for reduced to one scalar — was fitted and gating in
all three regions while `cell_scmd` was correctly declined in two: the same
quantity gated out under one name and admitted under another.

**Which columns are scMD-derived is read from the data, not from the name.** The
name `cell_composition_r2` is fixed by PI-locked configuration
(`config/aging.yml:axis.arms`, `config/gwas_negative_controls.yml`); the quantity
behind it is not. Module 04's MuSiC correction (also 2026-09-23) rebuilt it from
the RNA MuSiC proportions in every region and records the modality in
`cell_composition_r2_source`, keeping the scMD value beside it in
`cell_composition_r2_scmd`. So:

| Module 04 table | `cell_composition_r2` is | the arm gates in |
|---|---|---|
| up to `rra-AA-*-20260906` (accepted today) | scMD-derived | caudate only |
| after the MuSiC correction | MuSiC-derived | all three regions |

`age_functions.R:scmd_derived_feature_columns()` reads
`cell_composition_r2_source` and resolves the modality per column. A table with
no such column predates the correction and genuinely is scMD-derived, so the
constant is used as a fallback — announced in the run log and recorded as
`cell_composition_r2_source = unrecorded_assumed_dnam_scmd`, because it is an
assumption about data that cannot speak for itself. An **unrecognised** modality
is an error, never a silent "not scMD". The resolved modality is written to
`axis-tests.tsv`, `axis-arms-skipped.tsv`, `aging-decision.tsv` and the manifest,
so the audit trail names the modality that was actually adjusted for rather than
the one a constant assumed.

The judgement, since the column is a derived scalar and AGENTS.md §7.4 does ask
for "cell-composition-**associated** methylation properties": the adjective is
load-bearing. Where the proportions behind the R² do not track composition
(total-neuron ρ −0.03 DLPFC, −0.10 hippocampus, against 0.72 caudate) it measures
methylation variance shared with three PCs of something that is not composition,
and the only claim a gating boolean can license — "the gradient is not cell
composition" — is unavailable. §7.4 is asymmetric in the same way: MuSiC
adjustment is required in every region, scMD adjustment only "when the
integration gate passes", which `config/repeat_annotations.yml:334` states as
*caudate only, when the integration gate passes*. The selection is on modality;
it was never really about the name. The gate is still applied where the column is
**consumed**, not where it is built.

The **descriptive** annotation use of the column (stage 02's annotation table) is
unchanged and still fitted everywhere, flagged `scmd_derived_annotation` — which
now goes FALSE by itself once the column becomes MuSiC-derived. There the column
is the annotation being described, not an adjustment claiming to have removed
composition.

If the PI instead holds the column a methylation property, it should move out of
`same_n_members` for the reason `methylation_variance` is already non-gating —
a VMR's R² on donor PCs is partly a function of its total variance, which
contains the age variance — and the DLPFC reading changes the same way. The
readings only stay as sealed if the column is held to be both a methylation
property **and** a legitimate gate, in which case this code change is reverted.
`config/aging.yml` is `pi_locked` and was not edited; the arm's membership and
role are as locked.

The coverage gate `PASS_AGING_AXIS_COVERAGE` decides sealing. It is **not** a
success criterion: a null result seals.

### Cross-region reading (stage 05, `_m/combined/`)

The stage reads the Module 08 tiers. Caudate's tier is taken from
`config/region_donor_generalization.yml`.

- **Q1, region-general.** Module 09's H1 rule: supported in ≥ 2 regions,
  at least one of them not caudate. The endpoint is concordance, and no pooled
  p is formed. DLPFC and hippocampus share 115 of 118 donors, so they are not
  independent replicates.
- **Q2, identified difference (DLPFC vs hippocampus only).**
  - Variance = paired donor bootstrap over the donor union (within diagnosis)
    + a joint delete-one-chromosome jackknife of the difference.
  - "Context-dependent" is claimed only if the CI excludes zero in both
    `primary` and `cell_music`.
- **Q3, caudate.** Descriptive tier only. It is fitted and surfaced, and no
  caudate contrast is emitted (sequencing batch 3, AGENTS.md §8.1).

The decision token `aging_axis_association` takes one of these values:
- `REGION_GENERAL`
- `REGION_GENERAL_WITH_DLPFC_HIPPOCAMPUS_DIFFERENCE`
- `SINGLE_NONCAUDATE_REGION`
- `OPPOSITE_DIRECTION`
- `NOT_SUPPORTED`

Main-text vs supplement placement is a PI decision after the run.

## Catalog scope: a methylation PC tracks age

VMRs were selected on residual SD after regressing out methylation PCs 1–5.
Those PCs were computed **without** age adjustment (`01_vmr_catalog/_h/01_analyze.R`).
Stage 00 measures their correlation with age.

The smoke runs found a maximum |ρ| of **0.45** in DLPFC (chr11, PC5), **0.51**
in caudate (chr7, PC4) and **0.64** in hippocampus (chr17, PC4). Where a PC
tracks age, the catalog under-samples regions whose variability is mostly
age-driven. Results therefore describe age effects **within this catalog**,
which is the statement the manuscript may make. This is written down, not
corrected: correcting it would mean re-deriving the catalog.

## Interpretation constraints

These are carried on every decision row.

- The data are cross-sectional and postmortem. Write "age-associated
  methylation differences", **never** "methylation change with age". Age is
  confounded with birth cohort and survival.
- No causal, epigenetic-clock or environmental-determination claim. A null
  result is not evidence that low-control VMRs are environmentally determined
  (AGENTS.md §2.3), nor the reverse.
- The Module 02 score is never compared across regions. Only unit-free axis
  coefficients are.
- In DLPFC and hippocampus the only donor-level composition estimate is RNA
  MuSiC, because DNAm scMD fails its neuronal concordance gate there (ρ −0.03
  and −0.10).
  - Against the **currently accepted** Module 04 tables, `cell_composition_r2`
    is built from scMD proportions in all three regions, so since 2026-09-23 the
    arm that uses it is declined where the scMD gate fails and no longer gates
    the DLPFC or hippocampus reading. Against a **post-MuSiC-correction** table
    the column is MuSiC-derived, the arm is fitted in all three regions, and its
    estimate is a genuine composition adjustment there.
  - Residual composition confounding outside caudate is a DNAm-modality
    limitation either way: no DNAm-based composition estimate is usable there.
    The MuSiC-derived arm addresses RNA-estimated composition, not scMD's. That
    limitation therefore stands on its own and is **not** discharged by an arm
    whose covariate cannot measure the composition it names.

## Pipeline

| stage | script | output |
|---|---|---|
| 00 | `_h/00_new_run.R` | run, manifest, `diagnostic-methpc-age.tsv` |
| 01 | `_h/01_age_effects.R` | `vmr-age-effects.tsv`, `age-spec-summary.tsv`, `checkpoint/age-inputs.rds` |
| 02 | `_h/02_axis_test.R` | `axis-tests.tsv`, `axis-arms-skipped.tsv`, `axis-quartile-summary.tsv`, `annotation-age-associations.tsv`, bootstrap draws |
| 03 | `_h/03_apply_gates.R` | `gate-checks.tsv`, `gating-sensitivities.tsv`, `aging-decision.tsv` |
| 04 | `_h/04_finalize_run.R` | sealed run |
| 05 | `_h/05_cross_region_concordance.R` | `_m/combined/aging-*-AA.tsv`, including `aging-annotation-{associations,cross-region}-AA.tsv` |

To run:

    SMOKE_N=1 ./_h/submit_aging.sh AA dlpfc            # smoke: 50 draws, unlocked config
    ./_h/submit_aging.sh AA dlpfc                      # production, after the PI lock
    Rscript _h/05_cross_region_concordance.R --cohort AA   # after all three are accepted

Before any submission, run the smoke checks by hand:
`Rscript 09b_aging_application/tests/test_scmd_gate_consistency.R` (the scMD
gate, arm provenance, and the region reading reproduced from the sealed runs'
own `axis-tests.tsv`) and
`Rscript 09b_aging_application/tests/test_age_functions.R`. The latter has 19
checks, covering:
- matrix fits against `lm()`;
- donor alignment;
- stratified resampling;
- composite-null calibration and power;
- the mechanism check;
- determinism;
- cell PCs.

## Acceptance gate

`PASS_AGING_AXIS_COVERAGE` on all six coverage checks, with
`config/aging.yml` locked by the PI. Acceptance is a human entry below.

## Accepted runs

**The three runs below predate the 2026-09-23 scMD-gate correction and their
DLPFC reading is superseded.** `_m/` is immutable, so the correction takes effect
only in a new run. **What the rerun gives depends on which Module 04 table it
consumes, and the two answers differ**, so both are set out here.

### If Module 09b is rerun on the accepted `rra-AA-*-20260906` tables

`cell_composition_r2` is scMD-derived there, so its arm is declined outside
caudate and the reading follows from the sealed rows. Reading the sealed
`axis-tests.tsv` through the corrected stage-03 logic
(`tests/test_scmd_gate_consistency.R`) gives, with no refitting:

- **caudate: unchanged.** The scMD gate passes there, so the arm is legitimately
  fitted, and `cell_scmd` gates and fails at p 0.057 on its own.
  `PRIMARY_ONLY_FAILS_GATING_SENSITIVITY` stands.
- **DLPFC: `PRIMARY_ONLY_FAILS_GATING_SENSITIVITY` →
  `SUPPORTED_SURVIVES_GATING_SENSITIVITIES`.** `cell_composition_r2` was its only
  failing member; `cell_music` (−0.39, p 0.0096) and `controls_only` (72% of the
  primary) both survive.
- **hippocampus: unchanged, `NOT_SUPPORTED`.** The primary is not significant
  (p 0.28) and `cell_music` fails on its own.

Cross-region, that is one supported non-caudate region, not two, so the
association is still not region-general: the stage-05 token moves from
`NOT_SUPPORTED` to `SINGLE_NONCAUDATE_REGION`. The arithmetic above is a
projection from sealed numbers, not an accepted result; the accepted numbers
come from the rerun.

### If Module 04 is rerun first, and Module 09b consumes the MuSiC table

`cell_composition_r2` is then MuSiC-derived, the arm gates in **all three**
regions, and **none of the three readings can be projected**. The arm's estimate
is refitted on a different covariate: MuSiC PCs instead of scMD PCs, over more
donors (MuSiC covers 153/118/117 against scMD's 151/110/116, so DLPFC gains 8),
which can move it independently of the modality change. The sealed −0.17 (p 0.21)
for DLPFC is an scMD-based number and says nothing about the MuSiC-based one.

So DLPFC's reading is **not** guaranteed to flip in that state: it flips only if
the MuSiC-derived arm keeps the sign and reaches p < 0.05, which is an empirical
question for the rerun. The correction's guarantee is narrower and is the point
of it — whichever table is consumed, the arm gates exactly where its covariate
can support a composition claim, and `cell_composition_r2_source` on every
emitted row says which modality that was.

Run order therefore matters and is a PI decision: rerunning 09b alone tests the
axis against the currently accepted architecture, while rerunning Module 04 first
changes what the gating arm means. Module 09b does not require the Module 04
rerun — the gate is applied at consumption — so either order is valid, but the
two give different readings and the choice should be deliberate.

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| age-AA-caudate-20260919 | AA | caudate | vmrset-AA-caudate-937a41979978 | 2026-09-19 | Kynon J. Benjamin | PASS_AGING_AXIS_COVERAGE | Built at 0d69f432e on lgv-AA-caudate-rescore-20260913, rra-AA-caudate-20260906, tsc-AA-caudate-20260902; 153 donors, 11,251 VMRs, B = 1000. Primary −0.25 (p 3e-4), hypothesized direction; survives controls-only, MuSiC cell PCs, methylation variance, mappability, chromatin decomposition; fails the cell_composition_r2 arm (−0.07, p 0.19) and cell_scmd (p 0.057). Region reading PRIMARY_ONLY_FAILS_GATING_SENSITIVITY. Caudate is batch-confounded (descriptive tier). methPC–age max rho 0.51. |
| age-AA-dlpfc-20260919 | AA | dlpfc | vmrset-AA-dlpfc-856067dfe289 | 2026-09-19 | Kynon J. Benjamin | PASS_AGING_AXIS_COVERAGE | Built at 0d69f432e on lgv-AA-dlpfc-rescore-20260913, rra-AA-dlpfc-20260906, tsc-AA-dlpfc-20260902; 118 donors, 9,251 VMRs, B = 1000. Primary −0.33 (p 0.025), hypothesized direction; survives controls-only (72%), MuSiC cell PCs (−0.39, p 0.01), methylation variance, mappability, chromatin decomposition; fails the cell_composition_r2 arm (−0.17, p 0.21). Region reading PRIMARY_ONLY_FAILS_GATING_SENSITIVITY. methPC–age max rho 0.45. |
| age-AA-hippocampus-20260919 | AA | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-09-19 | Kynon J. Benjamin | PASS_AGING_AXIS_COVERAGE | Built at 0d69f432e on lgv-AA-hippocampus-rescore-20260913, rra-AA-hippocampus-20260906, tsc-AA-hippocampus-20260902; 117 donors, 9,166 VMRs, B = 1000. Primary −0.10 (p 0.28), hypothesized direction, not significant; controls-only −0.24 (p 0.02). Region reading NOT_SUPPORTED. methPC–age max rho 0.64, the largest of the three. |

Cross-region (stage 05, 2026-09-19): `aging_axis_association = NOT_SUPPORTED`.
Direction is concordant in all three regions and the top-quartile-control
VMRs carry 0.19–0.32 of the regional mean squared age effect, but the
`cell_composition_r2` gating arm removes the gradient everywhere. Reading:
the age-responsive low-control VMRs are the composition-sensitive ones.
Donor-level composition PCs in the age model do not remove it; the VMR-level
composition-sensitivity covariate does. Bulk data cannot separate an
age-related composition shift from a compartment that is cell-type-variable
and age-variable for the same reason; the manuscript may report the gradient
only with that qualifier. DLPFC–hippocampus difference −0.23, CI includes 0.

## Contract

`_h/` holds code, `_m/runs/{RUN_ID}/` holds immutable generated output, and
`tests/` is gitignored hand-run smoke checks. Configuration is shared at
`config/aging.yml`. See AGENTS.md §5.2 and §7.9.
