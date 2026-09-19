# 09b_aging_application — age-associated methylation along the local-genetic-control axis

**Status: implemented and smoke-tested 2026-09-19. `config/aging.yml` is
`pi_locked: false`; no production run has been submitted.**

## Question

Do VMRs with **weaker** local SNP control show **larger** age-associated
methylation differences?

This is an orthogonal, non-disease test of the pattern Module 09 found for
schizophrenia, where SCZ-linked VMRs sit lower on the local-genetic-control axis
in all three regions. If lower genetic control marks methylation that responds
more to aging, the axis carries biology beyond disease. The module is kept
lean on purpose: one axis test per region plus one cross-region stage. It is
not a second disease module.

The `09b` prefix follows the `01b` precedent. The module depends on 01, 02 and
04 and **not** on 09; it sits beside 09 as the second application of the axis.
It is not a sub-stage of 09.

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
| + `cell_composition_r2` | gating |
| + `methylation_variance` | non-gating |
| mappability ≥ 0.9 | non-gating |
| + H3K27me3, bivalent, H3K9me3 and quiescent fractions | non-gating decomposition: does the gradient ride on PRC2/bivalent age-hypermethylation or on Module 04's architecture? Attenuation links the two; it is not a failure. |
| signed β̂ (gain vs loss) | secondary, same inference; β̂ is unbiased, so no debiasing is needed |

### Region reading (stage 03)

The primary must be significant (p < 0.05) with a **negative** coefficient.
Then a strict conjunction applies:
- `cell_music`, `cell_scmd` and `cell_composition_r2` have the same n, so each
  must keep the sign **and** stay significant.
- `controls_only` keeps about 58% of donors, so it must keep the sign and at
  least 50% of the primary coefficient.

A gating member that cannot be fitted in a region is recorded as absent
evidence, not contrary evidence. `cell_scmd` outside caudate is the case in
point.

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
  - Module 04's `cell_composition_r2`, used by the gating arm, is built from
    scMD proportions in **all** regions, including those two.
  - Residual composition confounding therefore cannot be excluded outside
    caudate.

## Pipeline

| stage | script | output |
|---|---|---|
| 00 | `_h/00_new_run.R` | run, manifest, `diagnostic-methpc-age.tsv` |
| 01 | `_h/01_age_effects.R` | `vmr-age-effects.tsv`, `age-spec-summary.tsv`, `checkpoint/age-inputs.rds` |
| 02 | `_h/02_axis_test.R` | `axis-tests.tsv`, `axis-quartile-summary.tsv`, bootstrap draws |
| 03 | `_h/03_apply_gates.R` | `gate-checks.tsv`, `gating-sensitivities.tsv`, `aging-decision.tsv` |
| 04 | `_h/04_finalize_run.R` | sealed run |
| 05 | `_h/05_cross_region_concordance.R` | `_m/combined/aging-*-AA.tsv` |

To run:

    SMOKE_N=1 ./_h/submit_aging.sh AA dlpfc            # smoke: 50 draws, unlocked config
    ./_h/submit_aging.sh AA dlpfc                      # production, after the PI lock
    Rscript _h/05_cross_region_concordance.R --cohort AA   # after all three are accepted

Before any submission, run the smoke checks by hand:
`Rscript 09b_aging_application/tests/test_age_functions.R`. There are 19
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

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|

None. No production run has been submitted.

## Contract

`_h/` holds code, `_m/runs/{RUN_ID}/` holds immutable generated output, and
`tests/` is gitignored hand-run smoke checks. Configuration is shared at
`config/aging.yml`. See AGENTS.md §5.2 and §7.9.
