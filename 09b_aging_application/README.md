# 09b_aging_application — age-associated methylation along the local-genetic-control axis

**Status: accepted (AA, 2026-10-03).** `config/aging.yml` was locked by the PI
2026-09-19 and extended with the Fieller magnitude gate on 2026-10-03. The
accepted production runs are
`age-AA-{caudate,dlpfc,hippocampus}-20261003` (see **Accepted runs**, and
"### Superseded" for the `-20261001` and `-20260919` runs they replace, whose
DLPFC and caudate readings must not be quoted again). Cross-region token:
**`SINGLE_NONCAUDATE_REGION`** — the axis association is supported in DLPFC and
in DLPFC only. The 2026-09-19 `NOT_SUPPORTED` rested on a `cell_composition_r2`
built from a deconvolution invalid in two of three regions; Module 04's
2026-09-25 reacceptance made it MuSiC-derived and DLPFC survives gating on it.
Placement is settled: **supplementary figure, reported in the main text**
(PI, 2026-10-01; see "### Manuscript placement").

**The signed test and the magnitude are now licensed separately, and no region
carries both.** The 2026-10-03 round adds a Fieller denominator-stability gate
(T28 decision B, PI 2026-10-03): a *proportional* magnitude may be reported only
where the region's own mean debiased squared age effect -- the denominator the
outcome is divided by -- is itself separated from zero. DLPFC, the one supported
region, fails it at `den_z` 1.11, so its main-text statement is a **signed test
with an estimate and an interval, never a percentage**. Caudate (2.84) and
hippocampus (3.23) pass the gate and neither region is supported. The finding
itself is unchanged: DLPFC -0.328, p 0.0236, CI [-0.612, -0.044]. See
"## The 2026-10-03 magnitude-gate round".

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

The judgement, since the column is a derived scalar and the locked analysis plan does ask
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
  caudate contrast is emitted.

The decision token `aging_axis_association` takes one of these values:
- `REGION_GENERAL`
- `REGION_GENERAL_WITH_DLPFC_HIPPOCAMPUS_DIFFERENCE`
- `SINGLE_NONCAUDATE_REGION`
- `OPPOSITE_DIRECTION`
- `NOT_SUPPORTED`

Main-text vs supplement placement was a PI decision after the run, and was
made on 2026-10-01. See "### Manuscript placement".

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
  result is not evidence that low-control VMRs are environmentally determined, nor the reverse.
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

## The 2026-10-03 magnitude-gate round

`age-AA-{caudate,dlpfc,hippocampus}-20261003`, accepted 2026-10-03, replacing
`-20261001`. All four upstream run IDs and the `vmr_set_id` are identical in
every region, so nothing upstream moved: the round exists to install one gate.

**It adds no new primary number and moves no decision.** `estimate`, `n_vmrs`
and `n_chromosome_blocks` are bit-identical to `-20261001` in all three regions,
all three decision tokens are unchanged, and all three region readings are
unchanged. Only columns downstream of the donor bootstrap moved, which they must:
The locked analysis plan requires seeds derived from the run ID, so a new run ID necessarily
redraws them.

### Why the gate exists

Before this round the module had **no estimability check of any kind** on its
primary scale. `axis_estimate(scale = "relative_to_mean")` divided by
`mean(yy)` with no guard, so a region whose mean debiased squared age effect
happened to land at 1e-9 would have produced an enormous "proportional change"
and reported it. The proportional reading is only meaningful if the quantity it
is a proportion *of* is itself distinguishable from zero, and nothing was
checking that.

The gate is Fieller's denominator condition. With `den` the family mean and
`v22` its variance, a proportional magnitude is licensed only where
`den_z = den/sqrt(v22) > 1.96`. It gates the **magnitude and its interval only**
(PI, 2026-10-03). The signed test of whether the gradient is zero is unchanged
and stays on the ratio-scale combined SE, which T28 validated as calibrated.

### No region both supports the finding and permits a percentage

| region | region_reading | magnitude | den_z |
|---|---|---|---|
| caudate | PRIMARY_ONLY_FAILS_GATING_SENSITIVITY | REPORTABLE | 2.842 |
| DLPFC | SUPPORTED_SURVIVES_GATING_SENSITIVITIES | **NOT reportable** | 1.114 |
| hippocampus | NOT_SUPPORTED | REPORTABLE | 3.231 |

This is the substantive consequence of the gate for the manuscript. The one
region whose gradient survives gating is the one whose family mean is not
separated from zero, so it cannot carry a proportional magnitude; the two
regions whose means are well separated are the two where the gradient is not
supported. DLPFC's `primary_estimate_meaning` reads, in the sealed run:

> SIGNED TEST ONLY of whether the debiased squared age effect varies with
> local_snp_contribution_score_z; the magnitude is NOT reportable as a
> proportional change because the family mean is not separated from zero at
> alpha (family_mean_not_separated_from_zero_at_this_level, den_z 1.11)

**So 09b's main-text sentence is a signed statement about DLPFC with its
estimate and interval -- -0.328, p 0.0236, CI [-0.612, -0.044] -- and not a
percentage** (PI, 2026-10-03). Placement is untouched: the finding is reported
in the main text with its panel in the supplement, and only its scale changed.
Caudate's and hippocampus's rows keep the proportional wording, which they are
entitled to and which neither can use, because neither is supported.

### The gate reproduces its independent test fixtures

The gate was measured in a standalone test before the round ran, and the sealed
runs agree with it to within 4e-5 in all three regions -- the strongest
available evidence that the shipped code and the validated computation are the
same computation.

| region | test fixture | production emitted | diff |
|---|---|---|---|
| caudate | 2.8417 | 2.84168 | 1.7e-5 |
| DLPFC | 1.1142 | 1.11416 | 3.9e-5 |
| hippocampus | 3.2308 | 3.23077 | 2.7e-5 |

Its own error rates were measured rather than assumed, on 192 simulation cells
with a `no_effect` condition in which every per-VMR effect is exactly zero -- so
the true family mean is zero and no proportional effect exists:

- false "reportable" against a true zero denominator **0.0000** at n=152,
  0.0010 at n=48;
- max q95 `den_z` under a true zero denominator **≤ 1.03**, against the 1.96
  threshold;
- sensitivity **0.87-1.00** where the mean is well separated;
- `design` (per-VMR effects fixed vs redrawn) moves the gate not at all, which
  is correct -- the gate is a statement about the denominator, not about the
  effects.

### The pinned covariance is what keeps hippocampus reportable

`config/aging.yml` pins the gate's covariance as
`donor_jackknife_plus_chromosome_block`, and
`00_shared/axis_inference.R::require_magnitude_gate()` **refuses a run whose
config names anything else rather than defaulting to one**. The pin is
load-bearing, and in this module it is not hypothetical: the covariance choice
flips the gate across the 1.96 threshold in **5 of this module's 13 specs**.

| region | spec | bootstrap | pinned |
|---|---|---|---|
| caudate | cell_scmd \| base | 1.663 | 4.004 |
| hippocampus | primary \| base | 1.931 | 3.231 |
| hippocampus | primary \| cell_composition_r2 | 1.931 | 3.231 |
| hippocampus | cell_music \| base | 1.821 | 2.709 |
| hippocampus | controls_only \| base | 1.337 | 2.229 |

In every one of the five the donor-bootstrap covariance says NOT reportable and
the pinned covariance says reportable, including hippocampus's primary at 1.931
against a threshold of 1.96. Gating on the bootstrap would have imported the one
defect T28 established -- the donor bootstrap inflates the **absolute**-scale
variance, and `v22` is an absolute-scale variance -- into a decision, in 38% of
this module's specs. The pinned pairing is also the middle of three options, not
the most permissive: the donor jackknife alone is more liberal, so the
chromosome-block term is the conservative half of the pin.

### What T28 settled for this module, and the one caveat it leaves

T28 (signed 2026-10-03)
validated the inference this module uses. Four decisions were put and recorded:

- **D1 accepted.** The combined donor-bootstrap + delete-one-chromosome
  block-jackknife variance stays as the inference for the ratio scale. It is
  calibrated: type-I 0.027-0.046 against nominal 0.05 over 192 cells, and the
  best power of any calibrated variant. No config change.
- **D2 accepted.** The **raw**-scale rows are over-conservative and are not
  tests at 0.05. In this module that is the `signed_beta_age` secondary
  (`outcome_scale = raw`): its type-I error measured 0.000-0.003 against nominal
  0.05, so its p-values are too large and its intervals too wide. A rejection
  there is safe; a non-rejection means very little. No absolute-scale variance
  among those tested is calibrated, so an absolute-scale claim would need its
  own work.
- **D3 rejected.** The bootstrap draw guard admits a draw on the family mean
  being positive rather than separated from zero. The candidate replacement was
  measured and is worse -- power 0.041 against 0.246 at the realistic cell -- so
  the guard is unchanged and this is recorded as a known conservatism rather
  than fixed. The Fieller gate bounds the misuse of a near-zero denominator; it
  does not remove the inflated-draw effect on the gradient p.
- **D4 accepted.** `n_bootstrap` is not raised. For the quotient of a debiased
  statistic, raising B does not converge the series.

This module's primary p is **seed-stable at its own B**: 20 seeds at B = 1000
give 0.0227-0.0280, below alpha in 20 of 20, and 0.0243-0.0274 at B = 8000. The
`cell_composition_r2` gating arm behaves the same (0.0236-0.0290). The alpha
straddling T28 was raised to investigate is a Module 10 phenomenon and is not a
property of the variance method.

### Near-threshold rows in the annotation table

The re-runs surfaced one thing, and it is confined to
`annotation-age-associations.tsv` -- **not** `axis-tests.tsv`. Four DLPFC and
eight hippocampus rows sit close enough to q = 0.05 that the bootstrap seed
decides their call (q 0.044 → 0.052, 0.047 → 0.055; rows move by 0.004-0.016),
and fifteen and fourteen rows respectively lie in q [0.03, 0.08]. The counts at
q<0.05 therefore read 20/13/20 here against 20/17/26 in `-20261001`.

Two things this is not. It is **not a patch effect**: a patch-induced bias would
be directionally consistent across tables within a run, and it is not -- in
DLPFC the annotation p-values moved up in 75% of rows while the axis p-values
moved **down** in 78%. And it is **not 56 independent tests drifting**: those
rows share one set of bootstrap resamples, so they are correlated statistics,
and a single draw with slightly larger SEs moves most of them together. A sign
test on them assumes an independence that does not hold; the mean shift in p is
+0.002.

The locked analysis plan makes this table descriptive and explicitly outside the region
reading, so no decision depends on these rows. The consequence is a **writing
constraint**: do not report a near-threshold annotation as crisply significant
or crisply non-significant. Unlike the stage-B variance question, raising B
*would* genuinely fix this -- these are ordinary bootstrap SEs of a non-debiased
statistic, not a quotient of a debiased one -- and whether to spend that compute
on a descriptive table is open (PI, 2026-10-03: not decided).

### Known gap in this module's manifest

`manifest.tsv` does not record `relative_magnitude_gate` or its covariance,
while Module 10's does. The information **is** in the run, per-row in
`axis-tests.tsv` (`relative_magnitude_gate_reason`, `..._den_z`,
`..._gate_covariance`, `..._ci_lower/upper`,
`..._gate_n_donor_deletions_used`, `..._gate_decided_by_pi`), so a reader of the
results is fully served and a reader of the manifest alone is not. A code fix
for the next run, not a reason to re-run.

## Acceptance gate

`PASS_AGING_AXIS_COVERAGE` on all six coverage checks, with
`config/aging.yml` locked by the PI. Acceptance is a human entry below.

## Accepted runs

`00_shared/gates.R::read_accepted_runs()` takes the first `Accepted runs`
heading and stops at the next heading of any level, so the table must be the
first thing under this one. It was not: two `###` subsections sat above it, the
parser returned zero rows, and
`require_accepted_upstream("09b_aging_application", ...)` therefore refused every
consumer -- `11_integrated_manuscript_outputs/_h/11_supplementary_figures.R`
among them. Everything that was above the table is now below it.

The runs below postdate the 2026-09-23 scMD-gate correction. The three runs
they replace, and the one reading the rerun changed, are under
"### Superseded".

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| age-AA-caudate-20261003 | AA | caudate | vmrset-AA-caudate-937a41979978 | 2026-10-03 | Kynon J. M. Benjamin | PASS_AGING_AXIS_COVERAGE | Built at 9604a1529, git_dirty false, smoke_run FALSE, on vmrcat-AA-caudate-20260816, lgv-AA-caudate-rescore-20260913, rra-AA-caudate-20260925-a, tsc-AA-caudate-20260925-b -- all four upstreams and the vmr_set_id identical to the superseded `-20261001`, so the magnitude gate is the only thing that changed; 153 donors (88 control, 65 case), 11,251 VMRs modelled, 0 non-finite, 0 excluded, B = 1000 with 0 failures, 6/6 gate checks pass. Primary **-0.2517** (SE 0.0703, 95% CI -0.389 to -0.114, p 3.4e-04), hypothesized direction, on 11,204 VMRs over 22 chromosome blocks. `estimate`, `n_vmrs` and `n_chromosome_blocks` are bit-identical to `-20261001`; only bootstrap-derived columns moved. **Magnitude gate: REPORTABLE, den_z 2.842** (`fieller_denominator_stability`, covariance `donor_jackknife_plus_chromosome_block`) -- and not usable, because the region is not supported. Gating arms: controls_only -0.316 (p 2.6e-05, 126% of primary) PASS, cell_music -0.294 (p 5.7e-07) PASS, cell_composition_r2 -0.199 (p 9.8e-04) PASS, **cell_scmd -0.241 (p 0.0551) FAIL** -- the single failing member. Non-gating: methylation_variance -0.276 (p 4.1e-05), high_mappability -0.243 (p 1.6e-04, 8,196 VMRs), chromatin_decomposition -0.244 (p 4.9e-04), age_ge_25 -0.120 (p 0.199). Region reading **PRIMARY_ONLY_FAILS_GATING_SENSITIVITY**, region_supported FALSE. Quartile descriptive: top-control quartile carries 0.188 of the regional mean debiased squared age effect against 1.408 in the bottom quartile. 20 of 56 testable annotation associations at q<0.05 (the 57th row is untestable, fewer than 50 VMRs in the class; the superseded rows each counted one high). scMD integration gate PASS, so cell_scmd is legitimately fitted here and is the only region where it is. cell_composition_r2_source = rna_music. methPC-age max abs(rho) 0.513 at chr_7:PC4. Cross-sectional design; caudate is batch-confounded (descriptive tier) and no causal, clock or environmental-determination reading is licensed. |
| age-AA-dlpfc-20261003 | AA | dlpfc | vmrset-AA-dlpfc-856067dfe289 | 2026-10-03 | Kynon J. M. Benjamin | PASS_AGING_AXIS_COVERAGE | Built at 9604a1529, git_dirty false, smoke_run FALSE, on vmrcat-AA-dlpfc-20260816, lgv-AA-dlpfc-rescore-20260913, rra-AA-dlpfc-20260925-a, tsc-AA-dlpfc-20260925-b -- all four identical to the superseded `-20261001`; 118 donors (70 control, 48 case), 9,251 VMRs modelled, 0 non-finite, 0 excluded, B = 1000 with 0 failures, 6/6 gate checks pass. Primary **-0.3278** (SE 0.1449, 95% CI -0.612 to -0.044, p 0.0236; p was 0.0279 on the previous seed), hypothesized direction, on 9,214 VMRs over 22 chromosome blocks, bit-identical in estimate to `-20261001`. **Magnitude gate: NOT reportable, den_z 1.114** -- `family_mean_not_separated_from_zero_at_this_level`. `primary_estimate_meaning` therefore reads SIGNED TEST ONLY in the run itself, and **the main-text sentence for this module is a signed statement with the estimate and interval, not a percentage.** The finding is unaffected: the gate withholds the magnitude, not the result. Gating arms, **all fitted members survive**: cell_music -0.386 (p 0.0112) PASS, cell_composition_r2 -0.327 (p 0.0246) PASS, controls_only -0.236 (p 0.322 but 72.0% of the primary, against the locked `reduced_n_min_fraction: 0.5`; its rule is sign_and_min_fraction_of_primary, not significance, because the arm keeps only ~58% of donors) PASS; cell_scmd NOT FITTED, reason scmd_integration_gate_fails_in_region. Non-gating: methylation_variance -0.323 (p 0.0135), high_mappability -0.339 (p 0.0274, 5,702 VMRs), chromatin_decomposition -0.300 (p 0.0311), age_ge_25 -0.307 (p 0.143). Region reading **SUPPORTED_SURVIVES_GATING_SENSITIVITIES**, region_supported TRUE, unchanged from `-20261001`. Quartile descriptive: top quartile 0.237 against bottom 1.648. 13 of 56 testable annotation associations at q<0.05, against 17 in `-20261001`: four rows sat close enough to q = 0.05 that the bootstrap seed decides their call. That table is descriptive and outside the region reading, and the consequence is a writing constraint, not a decision -- see "### Near-threshold rows in the annotation table". scMD integration gate FAIL, cell_composition_r2_source = rna_music. methPC-age max abs(rho) 0.451 at chr_11:PC5. Cross-sectional design; "age-associated methylation differences", never "change with age". |
| age-AA-hippocampus-20261003 | AA | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-10-03 | Kynon J. M. Benjamin | PASS_AGING_AXIS_COVERAGE | Built at 9604a1529, git_dirty false, smoke_run FALSE, on vmrcat-AA-hippocampus-20260816, lgv-AA-hippocampus-rescore-20260913, rra-AA-hippocampus-20260925-a, tsc-AA-hippocampus-20260925-b -- all four identical to the superseded `-20261001`; 117 donors (69 control, 48 case), 9,166 VMRs modelled, 0 non-finite, 0 excluded, B = 1000 with 0 failures, 6/6 gate checks pass. Primary **-0.0977** (SE 0.0900, 95% CI -0.274 to 0.079, **p 0.278, not significant**), hypothesized direction, on 9,134 VMRs over 22 chromosome blocks, bit-identical in estimate to `-20261001`. **Magnitude gate: REPORTABLE, den_z 3.231** -- the strongest-separated denominator in the module, in the region whose gradient is not supported, and the clearest case that the gate and the finding are different questions. This is also where the pinned covariance is load-bearing: on the donor-bootstrap covariance the same row reads den_z 1.931 and would have been called NOT reportable. Gating arms: controls_only -0.243 (p 0.018, 249% of primary) PASS, cell_music -0.139 (p 0.165) FAIL, cell_composition_r2 -0.073 (p 0.414) FAIL; cell_scmd NOT FITTED, scmd_integration_gate_fails_in_region. Non-gating all null: methylation_variance -0.097 (p 0.284), high_mappability -0.017 (p 0.867, 5,776 VMRs), chromatin_decomposition -0.099 (p 0.258), age_ge_25 -0.049 (p 0.633). Region reading **NOT_SUPPORTED**, region_supported FALSE, unchanged. Quartile descriptive: top quartile 0.315 against bottom 1.145, the weakest separation of the three. 20 of 56 testable annotation associations at q<0.05, against 26 in `-20261001`: eight rows are seed-decided at q = 0.05, same writing constraint as DLPFC. methPC-age max abs(rho) **0.641** at chr_17:PC4, the largest of the three regions, so the catalog-scope caveat bites hardest here: a methylation PC removed before VMR calling tracks age, and this catalog therefore under-samples regions whose variability is mostly age-driven. That is a plausible contributor to the null and is not evidence against an age effect. |

**Cross-region (stage 05, rerun 2026-10-03 on the accepted `-20261003` runs):
`aging_axis_association = SINGLE_NONCAUDATE_REGION`**, with `citable = TRUE`,
`built_with_unaccepted_runs = FALSE`, `region_general = FALSE` and
`upstream_current TRUE` in all three regions. One region is supported (DLPFC),
one of them outside caudate, against the §7.9 rule's requirement of at least two
with at least one non-caudate. `_m/combined/` now holds this run. The token, the
three region readings and the identified-difference verdict are all unchanged
from the 2026-10-01 collation; what the tables gained is
`primary_relative_magnitude_reportable` and `primary_relative_magnitude_den_z`
per region, so a reader of `_m/combined/` meets the magnitude gate without
opening a run.

**The association is not region-general, and the module's qualifier is now about
count rather than composition — and, since 2026-10-03, about scale as well.** The superseded stage 05 returned
`NOT_SUPPORTED`, qualified by the `cell_composition_r2` arm removing the gradient
in every region. That qualifier does not survive (see "### Superseded"). What
replaces it is weaker and simpler: the gradient survives every fitted gating arm
in exactly one claim-eligible region, and the other claim-eligible region is
null. A single supported region cannot carry a region-general statement.

Q2, the identified difference, is **not claimed**:
`identified_difference_claimed = FALSE`. DLPFC minus hippocampus is **-0.230**
(95% CI -0.516 to 0.056) under `primary` and **-0.246** (CI -0.547 to 0.054)
under `cell_music`, unchanged to four decimals by the 2026-10-03 rerun; the rule needs the CI to exclude zero in both specs and it
excludes zero in neither. The difference is numerically unchanged from the
superseded -0.23, so the DLPFC flip moved the region reading without moving the
DLPFC-hippocampus contrast -- a reminder that a reading is a gating verdict, not
an effect size. The two regions share 115 of 118 donors (Jaccard 0.958) and are
not independent replicates; no pooled p is emitted
(`pooled_p_emitted = FALSE`).

Q3: caudate stays `descriptive_only` and no caudate contrast is emitted. `cross_region_raw_score_comparison_allowed = FALSE`,
`causal_interpretation_allowed = FALSE`,
`environmentally_determined_claim_allowed = FALSE`, and
`cross_sectional_design = TRUE` on the emitted row.

`manuscript_placement = supplementary_figure`, with
`manuscript_main_text_mention_required = TRUE` and
`manuscript_placement_provisional = TRUE`.

### Manuscript placement

**PI decision, 2026-10-01 (Kynon J. Benjamin), closing the open
"PI decision after the run" item.** 09b gets a **supplementary figure**, and the
finding is **reported in the main text**.

The reading is `SINGLE_NONCAUDATE_REGION`: one of three regions supported, no
region-general claim, and the DLPFC−hippocampus contrast does not exclude zero
(−0.230, CI −0.516 to +0.056). That does not carry a main figure.

Two things this decision is not:

- It is **not** permission to omit the result. Project-wide, every finding is
  stated in the main text even when its panel sits in the supplement;
  supplementary placement decides where the figure goes, never whether the
  result is reported. Module 11's `analysis-to-claim-matrix.tsv` is where that
  is enforced, so a 09b row must exist there with a main-text claim attached.
- It is **not** final in the upward direction. The placement may be revisited
  when the manuscript narrative is written, which cannot begin until every
  module has been rerun on current upstreams — Module 11 included. Revisiting
  placement does not reopen the result: the gating tokens are fixed by the
  accepted runs and a narrative does not move them.

Recorded machine-readably under `interpretation.manuscript_placement` in
`config/aging.yml`, and emitted on the `aging-cross-region-decision-AA.tsv` row
so the placement travels with the decision rather than living only in prose.

### Superseded

`age-AA-{caudate,dlpfc,hippocampus}-20261001` (accepted 2026-10-01, superseded
2026-10-03). Replaced by the magnitude-gate round and by nothing else. Their
upstreams, `vmr_set_id`, donor counts, VMR counts, primary estimates, decision
tokens and region readings are all identical to `-20261003`; what they lack is
the gate, so a reader of those runs could have quoted DLPFC's gradient as a
proportional change, which this project does not license. Their seed-dependent
columns differ by Monte Carlo error, and their annotation-table counts at
q<0.05 (20/17/26) differ from `-20261003`'s (20/13/20) for the reason in
"### Near-threshold rows in the annotation table". Nothing in them should be
re-quoted in preference to `-20261003`, and nothing in them is wrong.

`age-AA-{caudate,dlpfc,hippocampus}-20260919` (accepted 2026-09-19, superseded
2026-10-01). Superseded on two counts, and **the DLPFC and caudate readings
should not be quoted again.**

First, their upstreams moved. All three consumed `rra-AA-*-20260906` and
`tsc-AA-*-20260902`, both replaced by the 2026-09-25 acceptances. Modules 01 and
02 are unchanged, so the axis predictor is the same one. The locked analysis plan is not
retroactive: the 2026-09-19 acceptance was sound when it was made.

Second, and the reason a rerun was needed rather than a pointer, the
`cell_composition_r2` gating arm was **scMD-derived** in those runs. Outside
caudate the scMD integration gate FAILS, so in DLPFC and hippocampus that arm
was vetoing the reading on a covariate whose own modality is not licensed in
those regions. The 2026-09-23 stage-03 correction declines the arm where its
covariate cannot support a composition claim; with Module 04 reaccepted on
2026-09-25 the arm is MuSiC-derived and gates in all three regions.

**The composition qualifier does not survive, and this is the one place where
the module's conclusion changed rather than its accounting.** The superseded
cross-region reading was `aging_axis_association = NOT_SUPPORTED` with the
qualifier that "the age-responsive low-control VMRs are the
composition-sensitive ones" -- because `cell_composition_r2` removed the
gradient everywhere. On the MuSiC covariate it removes nothing in caudate
(-0.199, p 1.0e-03, against the superseded -0.07, p 0.19) and nothing in DLPFC
(-0.327, p 0.0289, against -0.17, p 0.21). The sealed -0.17 was an scMD-based
number and, as the projection below the old table said in advance, said nothing
about the MuSiC-based one.

So **DLPFC moves from `PRIMARY_ONLY_FAILS_GATING_SENSITIVITY` to
`SUPPORTED_SURVIVES_GATING_SENSITIVITIES`**, which was the flip the README
declined to project. Caudate stays `PRIMARY_ONLY_FAILS_GATING_SENSITIVITY` but
for a different and much narrower reason: its only failing member is now
`cell_scmd` at p 0.0562, in the one region where scMD legitimately gates.
Hippocampus is unchanged at `NOT_SUPPORTED`; its primary is p 0.269 and every
non-gating arm is null.

What did **not** change is worth recording. Direction is negative in all three
regions and in every arm of all three. Donor counts, VMR counts and the donor
checksums are identical to the superseded runs (153/118/117 donors;
11,251/9,251/9,166 VMRs), because Modules 01 and 02 did not move. The
`controls_only` arm survives in all three, and the top-control quartile still
carries far less of the regional mean squared age effect than the bottom
(0.188/0.237/0.315 against 1.408/1.648/1.145). The methPC-age diagnostics are
unchanged (0.513 / 0.451 / 0.641), so the catalog-scope caveat stands exactly as
written.

### Which Module 04 table was consumed, and why it mattered

Two projections used to sit here, one per run order, because the gating arm's
meaning depends on which Module 04 table supplies `cell_composition_r2`. The
question is now settled empirically and the projections are removed so they
cannot be quoted as results.

Module 04 was reaccepted first (`rra-AA-*-20260925-a`), so the arm is
**MuSiC-derived in all three regions** and `cell_composition_r2_source =
rna_music` on every emitted row. That is the branch the README said could not be
projected, for a reason worth keeping: the arm is refitted on a different
covariate over more donors (MuSiC covers 153/118/117 against scMD's 151/110/116,
so DLPFC gains 8), which can move it independently of the modality change. It
moved. The arm survives in caudate and DLPFC where the scMD-derived version
failed.

The correction's guarantee was always the narrower one, and it held: whichever
table is consumed, the arm gates exactly where its covariate can support a
composition claim, and `cell_composition_r2_source` on every row says which
modality that was.

## Contract

`_h/` holds code, `_m/runs/{RUN_ID}/` holds immutable generated output, and
`tests/` is gitignored hand-run smoke checks. Configuration is shared at
`config/aging.yml`. See the module README and its config.
