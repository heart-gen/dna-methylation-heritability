# 10_environmental_exploratory — exploratory environmental-factor analysis

**Status: three production runs sealed 2026-09-20, awaiting PI acceptance.
Supplement only, permanently.** `env-AA-{caudate,dlpfc,hippocampus}-20260920-a`
each returned `PASS_EXPLORATORY_COVERAGE` with 22/22 chromosomes reconciled and
zero computational failures. They supersede four earlier rounds; see
**Superseded runs**. Nothing here may be cited until the runs are entered in
**Accepted runs**.

Asks whether measured donor exposures associate with VMR methylation, and
whether that association varies along Module 02's continuous local-genetic-control
score. Both questions are exploratory and neither may enter the title, abstract,
primary groups, or main causal interpretation (AGENTS.md §2.3).

## Why this module exists, and what it is not

The v1 tree at `environmental-analysis/` is withdrawn — but for one specific
reason, not because the question was wrong. Every v1 endpoint read exposure
associations against `h2_category`, built in `00.prepare_data.R` from
`h2_unscaled` and `r_squared_cv`. Both metrics are retired (AGENTS.md §3), and
the genetically-anchored-vs-exposure-associated binary they encode is banned
outright (AGENTS.md §2.3). The scan survives; the classification it was read
against does not, and is not migrated.

PI decision 2026-09-10 closes the AGENTS.md §12 bullet "whether exposure
analyses remain supplemental or are removed": **retained, as an exploratory
supplemental module.**

**This module cannot show that a VMR is environmentally determined.** It has no
design that could. A null result here is consistent with an unmeasured exposure,
an exposure measured in the wrong window, or no exposure effect, and the module
is powered only for large effects of common exposures. `04_apply_gates.R` stamps
`environmentally_determined_claim_allowed = FALSE` on the decision row so this
travels with the result.

## The power limitation, made concrete

The demotion from primary is not rhetorical. Counted on the AA arm from
`inputs/phenotypes/_m/phenotypes-all.tsv`:

Two criteria, both prespecified: at least 20 donors in the minor class **and**
at most 15% missing. Minor-class count / fraction missing, on the AA arm.

Counted on the phenotype table after the age, region and race filters — the
*ceiling*. At run time the gate is applied to the Module 01 analysis set, which
is the intersection with the BSobj and the psam and is one donor smaller in
caudate (153) and one smaller in dlpfc and hippocampus (118, 117). The
eligibility outcome is the same either way; `results/exposure-eligibility.tsv`
carries the analysis-set numbers, which are the ones that governed a run.

| variable | caudate (n=154) | dlpfc (n=119) | hippocampus (n=118) | eligible |
|---|---|---|---|---|
| **smoking** (lifetime history) | 72 / 0.026 | 58 / 0.025 | 57 / 0.025 | all three |
| **nicotine** (toxicology) | 51 / 0.006 | 42 / 0.000 | 42 / 0.000 | all three |
| **education** | 41 / 0.110 | 30 / 0.134 | 31 / 0.136 | all three |
| **marital_status** | 37 / 0.104 | 27 / 0.118 | 28 / 0.119 | all three |
| **any_substance_nontobacco** | 26 / 0.006 | 19 / 0.000 | 19 / 0.000 | caudate only |
| antipsychotics | 46 / **0.234** | 37 / **0.269** | 37 / **0.271** | none |
| any_trauma_hx | 26 / **0.162** | 14 / **0.168** | 14 / **0.161** | none |
| ethanol | 13 / 0.013 | 9 / 0.000 | 9 / 0.008 | none |
| amphetamines | 7 / 0.013 | 5 / 0.017 | 5 / 0.017 | none |
| morphine | 5 / 0.006 | 4 / 0.000 | 4 / 0.000 | none |
| codeine | 4 / 0.006 | 4 / 0.000 | 4 / 0.000 | none |
| cocaine | 2 / 0.006 | 1 / 0.000 | 1 / 0.000 | none |
| fentanyl | 0 positives anywhere in the table | | | none |

Two exposures fail on **missingness, not on count**, and this is worth reading
carefully because it is not the reason the substance variables fail.
`antipsychotics` has 46 exposed donors in caudate but is unrecorded for 23% of
them; `any_trauma_hx` has 26 but is unrecorded for 16.2%, against a 15% ceiling.
Neither is a thin exposure — both are thinly *measured*. `any_trauma_hx` misses
the ceiling by 1.2 percentage points and is flagged `near_gate` in
`exposure-eligibility.tsv` for that reason. Relaxing the ceiling to admit them is
a PI decision and a config change, not something a run may do.

The practical consequence: the testable set is **tobacco, education and
marital_status in all three regions, plus any_substance_nontobacco in caudate**.
Substance use and trauma, the exposures the v1 analysis was largely about, do not
survive.

The gate lives in `config/environmental.yml` and
`01_build_exposure_matrix.R` writes every variable's counts and its pass/fail
reason to `results/exposure-eligibility.tsv`. v1 tested cocaine on two positive
donors in caudate and one in dlpfc; a "hit" there is a statement about one
person's methylome.

### Strata: annotation in this metadata is diagnosis-dependent

The pooled missingness above is not random, and treating it as random was the
error in the first version of this module. Annotation coverage differs sharply
by diagnosis, and it is the **control** stratum that is thin (AA caudate,
non-missing fraction):

| variable | Control (n=89) | Schizo (n=65) |
|---|---|---|
| antipsychotics | 61% | 98% |
| trauma history items | 76% | 88% |
| education | 83% | 97% |
| marital_status | 82% | 100% |
| smoking / nicotine | 96% / 100% | 100% / 98% |

Two exposures are affected in different ways.

**`antipsychotics` has zero exposed controls** (AA: Control 59 FALSE / 0 TRUE /
41 NA; Schizo 19 FALSE / 49 TRUE / 1 NA). Its entire contrast already lives
inside the case group, so a pooled test would present a within-case comparison
as a population-level exposure effect, routed through a collider. It is declared
`strata: [schizophrenia]` and is not testable pooled at any sample size.

**`lifetime_antipsych` is retired outright**, not for power: Control 76 FALSE /
0 TRUE, Schizo 0 FALSE / 68 TRUE. It *is* `primarydx`, recoded, and carries no
information beyond the diagnosis already in the locus model. It appears in
`config/covariates.yml:sensitivity_covariates.clinical`; it must never be used
as an exposure here.

**`any_trauma_hx` is a genuine unlock.** Pooled, it failed on 15.7% missing —
but that missingness is 23.6% in controls and 6.2% in cases. Inside the
schizophrenia stratum it clears the gate in caudate at 22 exposed of 65.
Stratifying fixes the *cause* of the failure rather than relaxing the limit.

So the module runs two strata: `all` (pooled, `primarydx` as covariate) and
`schizophrenia` (`primarydx` dropped — it is constant, and the collider path is
closed by the restriction itself rather than by adjustment). Each
**exposure × stratum pair is its own BH family**. A control stratum is
deliberately not declared: only `tobacco` clears the gate in controls, and
`tobacco` is already eligible pooled, so it would add a family without adding a
testable exposure.

#### The stratified floor of 15 is proportionally stricter, not looser

PI decision 2026-09-10. A count floor means different things at different sample
sizes, and 15 is the **higher** bar of the two as a fraction of the donors it
governs:

| gate | count | of | fraction |
|---|---|---|---|
| pooled | 20 | 153 (AA caudate) | 13.1% |
| stratified | 15 | 65 (Schizo caudate) | 23.1% |
| stratified | 15 | 48 (Schizo dlpfc/hippocampus) | 31.3% |

So a stratified exposure must clear a larger share of its stratum than a pooled
one clears of the full cell. That reasoning holds independently of any observed
count.

For the record, the number was set after the stratified counts were seen;
`config/environmental.yml` says so. The justification above does not rest on
them, and the pooled floor of 20 was locked before any count was seen. The floor
does not admit `antipsychotics` in dlpfc or hippocampus (11), `education` (12)
or `marital_status` (7) anywhere.

Eligible exposure × stratum pairs, AA caudate:

| stratum | eligible |
|---|---|
| `all` | smoking, nicotine, any_substance_nontobacco, education, marital_status |
| `schizophrenia` | nicotine (28), any_trauma_hx (22), smoking (19), antipsychotics (18) |

In dlpfc and hippocampus the schizophrenia stratum is n=48, where **`nicotine`
alone clears the gate** (22 of 48). That is a direct dividend of splitting it
from `smoking`: the union never cleared those cells, and `smoking` on its own
does not either (13 of 48).

### Drug-use composites

Regrouping the substance columns into pharmacologically coherent classes does
not recover power: opioids (codeine|morphine|fentanyl) reach 7/6/6, and
psychostimulants (cocaine|amphetamines) 9/6/6, across the three AA regions. Both
are recorded in `exposures.retired` with their counts so the attempt is not
silently repeated. Aggregation **across** mechanism is what crosses the gate, and
that yields two changes:

- **`smoking` and `nicotine` are deliberately NOT combined.** An earlier version
  unioned them as `tobacco`; that was wrong. Nicotine is toxicology — use at or
  near death — and smoking is lifetime history. Different exposure windows do
  not union. They are also empirically distinct: among AA donors they disagree
  for 28 of 164, with 25 smoking-history-positive but toxicology-negative (quit,
  or not recent) and 3 the reverse. Testing both is two questions, each with its
  own BH family.
- **`any_substance_nontobacco`** is a new eligible exposure in AA caudate only.
  It is a **burden indicator, not a pharmacological class** — an opioid-positive
  and an alcohol-positive donor share no mechanism — and no statement about any
  single drug is licensed by a result on it. dlpfc and hippocampus return 19
  against a gate of 20; `exposure-eligibility.tsv` flags that as `near_gate`
  rather than presenting a one-donor miss as a clean exclusion.

Resolved (cfc8dc145): `smoking` is lifetime history and `nicotine` is
toxicology at death. They disagree for 28 of 164 AA donors and are never
unioned; the earlier `tobacco` union is retired.

## Pipeline

| stage | script | what it does |
|---|---|---|
| 0 | `00_new_run.R` | opens the run; asserts the predictor, the FDR-family rule, the supplement-only flag and the eligibility keys before a directory exists; requires accepted 01, 02, 04 on one `vmr_set_id` |
| 1 | `01_build_exposure_matrix.R` | donor set from the accepted Module 01 run; derives the composites; applies the eligibility gate |
| 2 | `step_1_association.sh` → `02_vmr_exposure_association.R` | array 1-22: per-VMR `meth ~ exposure + covariates` |
| 3 | `02b_combine_associations.R` | reconciles the 22 tasks; **one BH family per exposure**, applied once |
| 4 | `03_control_axis_test.R` | the local-genetic-control axis models |
| 5 | `04_apply_gates.R` | coverage gate and interpretation flags |
| 6 | `05_finalize_run.R` | seals the run |
| 7 | `06_collate_regions.R` | collates the sealed runs into git-tracked tables under `_m/combined/`; **collates, never compares** |

    cd 10_environmental_exploratory/_m && mkdir -p logs
    ../_h/submit_environmental.sh AA caudate

Stages 0-6 run once per region. Stage 7 runs once, by hand, after all three are
sealed:

    Rscript 10_environmental_exploratory/_h/06_collate_regions.R --cohort AA

`SMOKE_N=1` permits unlocked keys and unaccepted upstreams; `DRY_RUN=1` prints
the job graph and submits nothing.

### The covariate model

`meth ~ exposure + age + sex + diagnosis`, built by
`00_shared/locus_io.R::load_locus_phenotype()` — the same reader
`load_observed_locus()` uses, so an exposure association conditions on exactly
what a local-genetic-control estimate conditions on (AGENTS.md §5.3).

No genotype PCs for the **arms**: population structure is removed upstream by
residualizing methylation on pooled snpPC1-3 during VMR discovery
(`config/covariates.yml:vmr_calling.genotype_pcs`). An estimation **cell**
carries `covs/genotype_pcs.tsv` from a within-group PCA and the reader picks it
up automatically. v1 used a global African-ancestry proportion here instead;
that is superseded.

### Stage B: the axis models

The v2 replacement for v1's `01.fishers_enrichment.py`. Three prespecified
models per eligible exposure:

1. **Primary, threshold-free:** `omega ~ score_z + vmr_length + cpg_count +
   cpg_density + gc_content + mappability + mean_methylation +
   methylation_variance`, where `omega` is the debiased exposure-explained sum of
   squares per donor, `(SS_exposure - df*sigma2_hat)/n`. Continuous on both
   sides, so the headline statement does not depend on where the FDR cut lands.
   `-log10(p)` was the primary until 2026-09-19 and is now a descriptive row
   only; see **The outcome scale couples to the score mechanically** for why.
2. Wilcoxon of `score_z`, FDR-significant vs not.
3. Logistic on the same indicator with the same covariates.

(2) and (3) are secondary — they need a threshold — and are skipped with a
recorded reason when either group holds fewer than 10 VMRs. Technical covariates
are joined from Module 04's `vmr-features.tsv` rather than recomputed, so GC
content and mappability are identical to the values Module 04 published.

## What a collaborator can read without Quest (`_m/combined/`)

`_m/runs/{RUN_ID}/` is immutable, gitignored and stays on Quest: the per-VMR
association table alone is 11 MB a region and the term matrix 20 MB. Stage 7
(`06_collate_regions.R`) writes the small tables that a manuscript actually
cites into `_m/combined/`, and those **are** tracked in Git, so nobody has to
re-run the SLURM chain to write up or check a number.

| file | rows | what it holds |
|---|---|---|
| `environmental-axis-reading-AA.tsv` | 19 | the stage B reading: one row per family with scale, `beta`, CI, `p`, `q`, the absolute sensitivity, `mean_omega_z`, the Fieller flag, the bootstrap inflation, the cell arm, and the Module 08 tier. A strict column subset of the table below — no number is derived here |
| `environmental-axis-per-region-AA.tsv` | 19 | every stage B column the runs emit (81) plus the collation's provenance columns, 90 in all, stacked |
| `environmental-vmr-associations-fdr-AA.tsv` | 12 | the stage A per-VMR rows that survive BH, with `vmr_id`, coordinates, `p_joint`, `fdr` and `omega`. All 12 are hippocampus (11 `nicotine`, 1 `education`); caudate and DLPFC have none |
| `environmental-fdr-families-AA.tsv` | 19 | per family: VMRs tested, VMRs significant, `min_p`, the BH alpha and method |
| `environmental-exposure-eligibility-AA.tsv` | 69 | every candidate exposure x stratum, its donor counts and class balance, and for the ineligible ones the reason. This is what explains why 19 families exist and not more |
| `environmental-decision-AA.tsv` | 3 | the per-region decision row, plus whether the Module 02 score it cites is still the accepted one |
| `environmental-gate-checks-AA.tsv` | 15 | each coverage check, observed against required |
| `environmental-collation-provenance-AA.tsv` | 3 | run ID, `sealed_at`, donor checksum, `vmr_set_id`, every upstream run ID, the run's git commit and the config SHA-256; plus the collation's own commit and the standing interpretation limits as text |

Three properties of these tables are worth stating because they are easy to
assume away:

- **No cross-region contrast.** `config/environmental.yml` sets
  `cross_region_comparison_allowed: false`, because caudate is sequencing batch 3
  and region is perfectly confounded with batch (AGENTS.md §8.1). Stage 7
  therefore stacks per-region rows and emits no pooled p, no region-general
  token and no between-region difference. Every row carries
  `cross_region_contrast_emitted = FALSE` and
  `regions_are_independent_replicates = FALSE`. The `tier` column marks caudate
  `descriptive_only` and the other two `claim_eligible`, as Module 08 requires.
- **FDR is not recomputed.** Every `q` is the one the sealed run wrote, from BH
  within that region's own family set. `fdr_recomputed_across_regions = FALSE`
  records it.
- **`citable` is a column, not a filename.** Module 09's stages 17/18 stamp
  `-UNACCEPTED` on their outputs because their *config* is unlocked;
  `environmental.yml` is locked, so that marker does not apply. What is missing
  is the PI's acceptance row. Until it exists, stage 7 refuses to run without
  `--allow-unaccepted-runs` and every table it writes carries `citable = FALSE`
  and `built_with_unaccepted_runs = TRUE`. Accepting the runs and re-running the
  stage overwrites the same paths, so the diff is exactly the acceptance flip.

## Interpretation constraints

- A null or negative axis result is **not** evidence that a VMR is
  environmentally determined (AGENTS.md §2.3).
- No heritable/non-heritable groups, no absolute PVE, no legacy metric.
  `forbidden_columns` is checked against every upstream table at runtime.
- **Region-specific only.** Caudate is sequencing batch 3 and region is
  perfectly confounded with batch (AGENTS.md §8.1), so no cross-region exposure
  contrast is emitted and none is licensed.
- `antipsychotics` clears the gate **only in caudate, and only inside the
  schizophrenia stratum** (18/46 split, 1.5% missing); in DLPFC and hippocampus
  the case stratum is n = 48 with a minor class of 11 against a floor of 15, and
  pooled it fails the 15% missingness ceiling in all three regions (22.9-26.5%).
  Where it is tested, the collider path is closed by the restriction rather than
  by adjustment — diagnosis is constant inside the stratum and is dropped from
  the locus model — and the row stays `collider_flagged` regardless.
- `any_substance_nontobacco` is a burden indicator across unrelated mechanisms.
- Supplement only, always. The decision row carries
  `main_text_retention = NEVER_SUPPLEMENT_ONLY`.

## First end-to-end run (smoke, 2026-09-10)

`env-smoke3-AA-caudate-20260910` ran all stages on the accepted upstreams. It is
a **smoke run and must not be cited**: it was driven stage-by-stage rather than
through the SLURM chain, so it does not satisfy AGENTS.md §9. It is recorded
because it establishes that the pipeline works, and because its shape is the
most informative thing the module has produced.

22/22 chromosomes reconciled; 11,341 VMRs per family; 9 exposure × stratum
families; `PASS_EXPLORATORY_COVERAGE`.

**Stage A is null in every family.** Zero FDR-significant VMRs, all nine. The
smallest p anywhere was 4.7 × 10⁻⁶ (nicotine, pooled) across 11,341 tests —
nothing close to surviving BH.

**Stage B: every pooled gradient vanishes inside the case stratum.**

| exposure | stratum | beta | p |
|---|---|---|---|
| nicotine | all | −0.0368 | 5.2 × 10⁻¹⁰ |
| smoking | all | −0.0325 | 2.5 × 10⁻⁸ |
| marital_status | all | 0.0181 | 1.9 × 10⁻⁴ |
| any_substance_nontobacco | all | 0.0164 | 8.1 × 10⁻⁴ |
| education | all | 0.0088 | 0.079 |
| smoking | schizophrenia | 0.0021 | **0.67** |
| nicotine | schizophrenia | −0.0018 | **0.72** |
| antipsychotics | schizophrenia | 0.0052 | **0.29** |
| any_trauma_hx | schizophrenia | 0.0039 | **0.43** |

**This is a negative control that worked, not a finding.** Three observations
line up:

1. Stage A found nothing, so Stage B is describing the shape of a p-value
   distribution that is entirely noise.
2. The exposures with the strongest pooled gradients are the
   **diagnosis-correlated** ones — smoking (72% of cases vs 28% of controls) and
   nicotine (46% vs 25%) — while `education`, the least diagnosis-linked, is the
   weakest at p = 0.079.
3. Every gradient collapses to null when diagnosis is held constant, with
   coefficients falling by more than an order of magnitude.

The straightforward reading is that the pooled gradient tracks *diagnosis*, or
something technical correlated with it, rather than exposure. It is not evidence
that exposure effects concentrate in low-genetic-control VMRs, and it must not
be written up as such — that is precisely the inference AGENTS.md §2.3 forbids.
The competing mundane explanation, n=65 against n=153, is also not excluded by
this design.

The grouped tests could not run in any family
(`fewer_than_10_vmrs_in_a_group`), so no row above has a second line of
evidence.

## First production runs (2026-09-19) — superseded, kept for the record

**These numbers are on the retired `-log10 p` scale and must not be cited.** The
section is kept because the two design changes below were found in these runs,
and because the contrast between the two scales is the evidence for the change.


`env-AA-{caudate,dlpfc,hippocampus}-20260919`, submitted through the SLURM chain
on the accepted `lgv-AA-*-rescore-20260913` score and `rra-AA-*-20260906`
features. All three: 22/22 chromosomes expected / completed, 0 excluded, 0
QC-failed, 0 failed, 0 unaccounted; all five coverage checks pass;
`PASS_EXPLORATORY_COVERAGE`.

| region | donors | VMRs per family | eligible exposure×stratum | stage A pairs at q<0.05 | stage B families at q<0.05 |
|---|---|---|---|---|---|
| caudate | 153 | 11,251 | 9 (4 in cases) | **0** | 4 of 9 |
| dlpfc | 118 | 9,251 | 5 (1 in cases) | **0** | 3 of 5 |
| hippocampus | 117 | 9,166 | 5 (1 in cases) | **12** | 4 of 5 |

**Stage A.** Null in caudate and DLPFC, as in the smoke run. Hippocampus returns
12 FDR-significant VMR × exposure pairs, 11 of them `nicotine` (3 pooled, 8
inside the case stratum at n = 48) and one `education`. `chr8:47762980-47763627`
is the only VMR significant in both strata (q = 0.010 pooled, 0.003 in cases).
Twelve pairs out of five families of ~9,200 tests is a thin result in the one
region with the smallest exposed case count, and it is a supplemental
observation, not a finding.

**Stage B, the axis.** Only the threshold-free primary reports: the Wilcoxon and
logistic forms were skipped in **every** family in **every** region
(`fewer_than_10_vmrs_in_a_group`), which is why the primary was specified
threshold-free.

| exposure | stratum | caudate β (q) | dlpfc β (q) | hippocampus β (q) |
|---|---|---|---|---|
| nicotine | all | −0.036 (1.1e-8) | −0.026 (1.5e-4) | −0.052 (1.2e-12) |
| smoking | all | −0.033 (6.3e-8) | −0.002 (0.66) | −0.017 (0.0078) |
| marital_status | all | +0.017 (6.1e-4) | −0.022 (1.5e-4) | −0.006 (0.24) |
| education | all | +0.008 (0.098) | +0.009 (0.096) | **+0.023 (1.8e-5)** |
| any_substance_nontobacco | all | +0.016 (0.0014) | not eligible | not eligible |
| nicotine | schizophrenia | −0.0005 (0.93) | **−0.042 (6.5e-9)** | **−0.088 (8.2e-30)** |
| smoking | schizophrenia | +0.003 (0.77) | not eligible | not eligible |
| antipsychotics | schizophrenia | +0.004 (0.74) | not eligible | not eligible |
| any_trauma_hx | schizophrenia | +0.005 (0.74) | not eligible | not eligible |

β is the change in the per-VMR exposure −log₁₀p per SD of
`local_snp_contribution_score_z`, so these are shifts of a few hundredths of a
−log₁₀p unit across the whole axis. They are statistically sharp because there
are ~9,000-11,000 VMRs per family and numerically negligible.

**Two things changed relative to the caudate smoke run.** Caudate reproduces it
exactly — every pooled gradient collapses inside the case stratum. DLPFC and
hippocampus do **not**: `nicotine` in cases is the *strongest* gradient in both
(the hippocampal coefficient is 1.7× its own pooled value), so "the pooled
gradient is diagnosis mixing" is a caudate statement and does not generalize. And
`education`, the least diagnosis-linked exposure, is FDR-significant in
hippocampus with the **opposite** sign to nicotine, so the axis gradients do not
share a direction across exposures.

### The outcome scale couples to the score mechanically

This must be settled before any of stage B is written up. The primary outcome is
−log₁₀p of a per-VMR exposure model that contains **no SNP term**, so a
high-control VMR carries its local genetic variance in that model's residual.
That inflates its SE and deflates its −log₁₀p, which produces a negative
coefficient on `score_z` with no exposure biology involved. `09b_aging_application`
hit the identical problem on 2026-09-19 and rejected every SE-dependent outcome
(|β|, rank|β|, |t|, −log₁₀p) in favour of the debiased `β̂² − SE²`, with donor
bootstrap plus block-jackknife variance (AGENTS.md §7.9).

#### What the diagnostic found (2026-09-19, hand-run, `tests/debiased_axis_firstlook/`)

`vmr-exposure-association-terms.tsv` already carries per-VMR `beta` and `se`, so
for the five **1-df** exposures the debiased outcome `ω = β̂² − SE²` is computable
from the sealed runs with no refit. The two 2-df exposures (`education`,
`marital_status`) are not: they need the anova sum of squares and σ̂², which
stage 2 does not emit. Three things came out of it, and they reorder the problem.

**1. The gradient is not an artifact of the p-value scale.** Every family that is
FDR-significant on `−log₁₀p` is significant on `ω` with the same sign, and the two
outcomes order the VMRs almost identically (Spearman 0.91-0.96). On the
interpretable scale the effect is **large**: the mean debiased squared exposure
effect falls by 23-49% per SD of score (caudate nicotine −46%, smoking −47%;
hippocampus nicotine −38%, within cases −49%; DLPFC nicotine −23%, within cases
−27%). The "numerically negligible" reading of the `−log₁₀p` coefficients was an
artifact of that scale, not a property of the association.

**2. The residual-variance channel is real and carries much of the pooled
result.** Within a family the design is shared, so `SE²` *is* the exposure model's
residual variance up to a constant. It rises with `score_z` with total
methylation variance held fixed, in every pooled family (p 1e-13 to 1e-28) — the
local genetic variance is demonstrably sitting in that residual. Holding it fixed
attenuates the `ω` gradient by 29-57% in the pooled families and removes it in
two (`nicotine@all` in DLPFC, `smoking@all` in hippocampus). It does **not** touch
the strongest result: `nicotine` within cases attenuates 3% in hippocampus
(p = 4e-29) and gets *stronger* in DLPFC. This decomposition is indicative only —
`SE²` is a component of `ω`, so conditioning on it partly conditions on the
subtracted term. A clean version needs σ̂² from the covariate-only null model,
which is the same extra column the 2-df exposures need.

**3. Stage B's inference is anti-conservative, independently of the outcome.**
`03_control_axis_test.R` reports plain `lm()` p-values over 9,000-11,000 VMRs
treated as independent observations. A delete-one-chromosome weighted block
jackknife — the construction Modules 06, 08 and 09b use — gives SEs a **median
1.28× larger** (up to 2.0×), moving p-values by two to five orders of magnitude.
Only one verdict at 0.05 changes among the significant families (caudate
`any_substance_nontobacco` on `ω`, 0.012 → 0.069), so the direction of the result
holds; but **the q-values now in the run's tables are optimistic**, and this
applies to what is already reported rather than to a proposed replacement. The
jackknife is only the VMR-level half — the donor-level half, which every VMR
shares through the same 117-153 donors, is not in it at all.

#### What followed: PI approval 2026-09-19

Both changes were approved by the PI on 2026-09-19 and are implemented.

- **Primary outcome** is now `debiased_partial_ss`, the debiased
  exposure-explained sum of squares per donor, `(SS_exposure − df·σ̂²)/n`,
  computed in stage 2 where the model is fitted. It is unbiased for the exposure
  term's contribution whatever the SE, zero in expectation under the null, and
  negative values are retained. For a 1-df exposure it equals `β̂² − SE²` up to a
  factor constant within a family; the generalized form is what lets the two 2-df
  categorical exposures (`education`, `marital_status`) use the same endpoint.
  `−log₁₀p` is retained as a **descriptive** row with its own BH family and
  `mechanically_biased_toward_hypothesis = TRUE`.
- **Inference** is donor-bootstrap variance + delete-one-chromosome weighted
  block-jackknife variance, B = 2000, seeded per family from the run ID. The
  generic machinery moved out of 09b into `00_shared/axis_inference.R`
  (`resample_rows`, `prepare_axis`, `axis_estimate`, `block_jackknife_se`,
  `combined_inference`, plus `fit_scalar_matrix`/`fit_term_matrix` and the two
  debiasing helpers); 09b sources it and
  `09b_aging_application/tests/test_shared_axis_equivalence.R` checks the shared
  versions against the definitions 09b carried at `abb0c6789`, so no accepted 09b
  number moved. The OLS p-value is still emitted, as
  `primary_p_ols_vmrs_independent`, so the difference stays visible.
- Stage 2 now writes a donor × VMR phenotype checkpoint, because the bootstrap
  refits every VMR 2000 times and cannot re-read one phenotype file per draw.
  `_h/exposure_model.R` holds the stratum filter, the covariate list and the
  design matrix **once**, called by both stages, so the bootstrap cannot drift
  from the model whose estimate it is putting a SE on. Stage 3 asserts that its
  matrix refit reproduces stage 2's per-VMR `ω` and stops if it does not.
- `00_new_run.R` refuses a config that names any other primary outcome or
  inference, rather than defaulting, the way
  `gates.R::donor_group_inference_policy()` refuses a widened donor-group policy.

### On what scale the gradient is reported (PI 2026-09-20)

`ω` is an exposure-explained variance, and its *level* is a nuisance: it varies
about 5× between chromosomes and about 6× between donor draws. Dividing by the
family mean before fitting makes the coefficient a **proportional** change per SD
of score, and that is the estimand the earlier runs reported. The question the PI
put was whether the relative scale should be kept at all, and if so how to decide
when its denominator is too close to zero to divide by. The old answer was a
hand-set `relative_scale_min_mean_z: 2` guard; it was wrong twice over, because
its z treated 9,000-11,000 correlated VMRs as independent and it admitted the very
family it was written to exclude. It is **retired**, and `00_new_run.R` hard-stops
on a config that still names it.

The design now:

- **The primary is the proportional gradient**, estimated as a ratio functional:
  every bootstrap draw and every deleted chromosome recomputes `β/mean(ω)` on its
  own rows, which is the standard variance treatment for a ratio estimator.
- **The absolute gradient rides on every row** as a sensitivity —
  `absolute_beta`, `absolute_se`, `absolute_p`, `absolute_role =
  sensitivity_not_gating`. It gates nothing. A reader can see both scales rather
  than take the module's word for one.
- **A Fieller interval for the percentage is reported and never gates**
  (`fieller_bounded`, `fieller_role = informational_not_gating`). On these data it
  is unbounded in every family: the family mean is separated from zero by only
  `mean_omega_z` 0.07-1.60, so no percentage is formally identifiable. That is a
  statement about the denominator's own weakness, not about the primary.
- **No threshold is set anywhere.** A gate on the denominator's separation from
  zero excludes every family, including all three findings; and a gate on the
  ratio's own interval (`CI width > 2|β̂|`) is algebraically `p > 0.05`, so it
  would only restate the p-value already on the row. The one boundary left is
  definitional: `mean(ω) ≤ 0` means no detectable exposure contribution, so the
  ratio has no sign and that family reports the absolute gradient
  (`primary_scale = raw`).
- **`bootstrap_mean_omega_inflation` is on every relative-scale row.** Resampling
  donors with replacement duplicates donors, which inflates an exposure-explained
  sum of squares: the factor runs 2.0×-38.5× here. Where it exceeds 1 the ratio's
  **bootstrap** variance is likely optimistic, so these p-values may be too small.
  That is not resolved; it is made auditable. The same bootstrap construction is
  in 09b's accepted runs, so it reaches there too. Settling it needs a null and
  positive-control simulation that has not been run.

**A defect found on 2026-09-20 and the reason `env-AA-*-20260920` is superseded.**
`block_jackknife_cov()` computed its pseudo-values as `hj * full`, recycling a
length-`n_blocks` vector against a length-2 one instead of forming the outer
product. With 22 chromosomes — an even multiple of 2 — R issues **no warning**, so
it was silent. It corrupted only the joint covariance, hence the `absolute_*`,
`fieller_*` and `mean_omega_*` columns; the primary goes through
`block_jackknife_se()`, a scalar path with no recycling, so no primary estimate,
SE, p or q was ever affected. It was caught because on a `raw`-scale family the
primary *is* the absolute estimand, yet their jackknife SEs disagreed ~8×. The
invariant now has a regression test,
`tests/test_jackknife_cov_agrees.R`, which fails on the old code at 22 blocks.
Note that `code/` snapshots `_h/` and `config/` but **not** `00_shared/`, which is
how a shared-code defect reached three sealed runs with no provenance trail.

### Limitation for the discussion: the variance budget

This is not a task. No analysis in this module can remove it, and none is
planned. It is the paragraph to carry into PI summaries and into the manuscript
Discussion wherever the Stage B gradient is mentioned:

> The axis contrast holds each VMR's total methylation variance fixed, so a VMR
> higher on the local-genetic-control axis has, by construction, less non-genetic
> variance left for any exposure to move. A negative gradient is therefore close
> to arithmetically expected wherever a real exposure effect exists, and it is
> not evidence that exposure effects are concentrated among VMRs with weaker
> local genetic control. Separating the two readings would require expressing
> each exposure effect relative to that VMR's non-genetic variance, and that
> denominator is an absolute local-PVE estimate — a quantity this project
> retired, because its simulation calibration passed for relative ordering and
> failed for absolute PVE. Debiasing the outcome removes the statistical half of
> the problem, the SE inflation that makes a high-control VMR look less
> exposure-responsive than it is; it cannot remove the arithmetic half. The
> gradient is therefore reported in absolute units, as a description of where
> exposure-associated variance sits along the axis, and never as evidence that a
> VMR is environmentally determined.

Recorded machine-readably as `interpretation.variance_budget_limitation` in
`config/environmental.yml`, and in AGENTS.md §7.10. It applies equally to §7.9's
aging axis, which asks the same question with age in place of exposure.

**One open PI item follows from it.** `methylation_variance` is in this module's
prespecified `axis_covariates`, and holding it fixed is exactly what activates the
budget arithmetic. 09b faced the same choice and **excluded** total variance from
its primary, on the ground that the outcome is part of it. The covariate set is
prespecified and `pi_locked`, so it was left as locked here rather than changed
alongside the outcome; whether Module 10 should carry a no-`methylation_variance`
arm is a PI decision, and it is cheap to add now that the machinery exists.

## Current production runs (`env-AA-*-20260920-a`, sealed 2026-09-20)

On the accepted `lgv-AA-*-rescore-20260913` score and `rra-AA-*-20260906`
features, unchanged from the 2026-09-19 round. All three: 22/22 chromosomes
expected / completed, 0 excluded, 0 QC-failed, 0 failed, 0 unaccounted; all five
coverage checks pass; `PASS_EXPLORATORY_COVERAGE`.

| region | run_id | donors | VMRs (tested) | eligible exposure×stratum | stage A pairs at q<0.05 | stage B families at q<0.05 |
|---|---|---|---|---|---|---|
| caudate | `env-AA-caudate-20260920-a` | 153 | 11,251 (11,204) | 9 (4 in cases) | **0** | 1 of 9 |
| dlpfc | `env-AA-dlpfc-20260920-a` | 118 | 9,251 (9,214) | 5 (1 in cases) | **0** | 0 of 5 |
| hippocampus | `env-AA-hippocampus-20260920-a` | 117 | 9,166 (9,134) | 5 (1 in cases) | **12** | 2 of 5 |

**Stage A is unchanged** from 2026-09-19, and necessarily so: it carries no
bootstrap, so it is deterministic given the same inputs. Null in caudate and
DLPFC. Hippocampus returns the same 12 FDR-significant VMR × exposure pairs, 11
`nicotine` and one `education`, with `chr8:47762980-47763627` again the only VMR
significant in both strata (q = 0.0100 pooled, 0.0027 in cases). Twelve pairs out
of five families of ~9,200 tests, in the region with the smallest exposed case
count, remains a supplemental observation and not a finding.

**Stage B.** The Wilcoxon and logistic forms were again skipped in every family in
every region (`fewer_than_10_vmrs_in_a_group`), so only the threshold-free primary
reports. `β` is the **proportional** change in the debiased exposure-explained sum
of squares per SD of `score_z`, as a fraction of the family mean — multiply by 100
for a percentage. Families whose mean `ω` is at or below zero cannot form that
ratio and report the absolute gradient instead (`primary_scale = raw`), marked ᴿ.

| region | exposure | stratum | primary β | p | q | absolute p | mean_omega_z | boot infl |
|---|---|---|---|---|---|---|---|---|
| hippocampus | nicotine | schizophrenia | **−0.491** | 1.1e-5 | **1.1e-5** | 0.081 | 1.60 | 2.0× |
| hippocampus | nicotine | all | **−0.376** | 3.3e-4 | **0.0013** | 0.21 | 1.49 | 2.6× |
| caudate | nicotine | all | **−0.462** | 0.0017 | **0.0085** | 0.41 | 0.78 | 5.1× |
| hippocampus | smoking | all | −0.322 | 0.029 | 0.057 | 0.53 | 0.86 | 4.1× |
| dlpfc | nicotine | schizophrenia | −0.267 | 0.054 | 0.054 | 0.45 | 0.76 | 3.2× |
| caudate | smoking | all | −0.474 | 0.030 | 0.074 | 0.43 | 0.70 | 5.8× |
| dlpfc | nicotine | all | −0.232 | 0.134 | 0.269 | 0.60 | 0.66 | 4.0× |
| dlpfc | marital_status | all | −2.411 | 0.106 | 0.269 | 0.58 | 0.07 | 38.5× |
| dlpfc | smoking | all | −0.113 | 0.900 | 0.900 | 0.97 | 0.09 | 27.0× |
| hippocampus | education | all | ᴿ +1.06e-5 | 0.242 | 0.323 | 0.242 | −0.49 | — |
| hippocampus | marital_status | all | ᴿ −4.98e-6 | 0.695 | 0.695 | 0.695 | −0.04 | — |
| caudate | marital_status | all | ᴿ +5.85e-6 | 0.345 | 0.542 | 0.345 | −0.41 | — |
| caudate | any_substance_nontobacco | all | ᴿ +2.16e-6 | 0.533 | 0.542 | 0.533 | −0.27 | — |
| caudate | education | all | ᴿ +3.18e-6 | 0.542 | 0.542 | 0.542 | −0.34 | — |
| dlpfc | education | all | ᴿ +4.23e-6 | 0.763 | 0.900 | 0.763 | −0.58 | — |
| caudate | smoking | schizophrenia | ᴿ +1.63e-6 | 0.862 | 0.980 | 0.862 | −0.09 | — |
| caudate | any_trauma_hx | schizophrenia | ᴿ −1.02e-6 | 0.928 | 0.980 | 0.928 | −0.06 | — |
| caudate | antipsychotics | schizophrenia | ᴿ +2.39e-7 | 0.979 | 0.980 | 0.979 | −0.26 | — |
| caudate | nicotine | schizophrenia | ᴿ +2.49e-7 | 0.980 | 0.980 | 0.980 | −0.14 | — |

**Three families survive FDR, all `nicotine`, and all negative.** Hippocampus in
the case stratum (−49%, q = 1.1e-5) and pooled (−38%, q = 0.0013), and caudate
pooled (−46%, q = 0.0085). The `cell_composition_r2` arm barely moves any of them
— hippocampus cases −0.483 (p = 1.9e-6), hippocampus pooled −0.370 (p = 1.4e-4),
caudate pooled −0.440 (p = 0.0040) — so the gradient is not cell composition
restated. Three more families are marginal and resolve nothing: hippocampus
`smoking` (q = 0.057), dlpfc `nicotine` in cases (q = 0.054) and caudate `smoking`
(q = 0.074).

**What the sensitivity columns say, and they are not decoration.** `absolute_p` is
above 0.05 in **every** family in all three regions (`n_absolute_p_below_alpha =
0`), and the Fieller interval is unbounded in every family
(`n_fieller_bounded = 0`), because `mean_omega_z` never reaches 1.96 — its maximum
anywhere is 1.60. So the *level* of exposure-explained variance is not resolved by
these data at any locus, and **no percentage in the table above is formally
identifiable as a percentage**; the ratio is what replicates, not the scale it is
a ratio of. On `raw`-scale families the primary and the absolute estimand are the
same quantity, and their p-values agree exactly, which is the check that the two
paths compute one thing.

**Four things a reader must not do with this table.**

1. Do not read `−0.491` as "49% less exposure-explained variance" without the
   caveat above. The point estimate is a ratio; its denominator is not separated
   from zero.
2. Do not quote a percentage from a `raw` row at all. Its family mean is at or
   below zero, so there is no percentage to quote.
3. Do not read `arm_attenuation` where the primary is null. It is a ratio of two
   coefficients, so it explodes as the denominator approaches zero: in these runs
   it is −8.6 for caudate `antipsychotics` (p = 0.98), +4.2 for caudate
   `any_trauma_hx` (p = 0.93), +3.3 for dlpfc `smoking` (p = 0.90) and +2.2 for
   caudate `nicotine` in cases (p = 0.98). None of those means anything. In the
   three families that survive FDR it is 0.016-0.048, i.e. the arm changes almost
   nothing, which is the only place the column is worth reading.
4. Do not treat the negative sign as evidence that exposure effects concentrate
   at weakly controlled VMRs. **Limitation for the discussion: the variance
   budget** explains why a negative gradient is close to arithmetically expected,
   and that limitation is permanent.

## Superseded runs

| run_id | sealed | superseded because |
|---|---|---|
| `env-AA-{region}-20260919` | 2026-09-19 | primary was `-log10 p`, which couples to the score mechanically, and inference was an OLS p-value over 9,000-11,000 correlated VMRs |
| `env-AA-{region}-20260919-a` | 2026-09-19 | first debiased-outcome round; the `relative_scale_min_mean_z` guard it shipped with was invalid and inert |
| `env-AA-{region}-20260919-b` | 2026-09-19 | the guard fix did not fire; reproduces `-a` exactly up to bootstrap noise |
| `env-AA-{region}-20260920` | 2026-09-20 | `block_jackknife_cov()` recycling defect corrupted every `absolute_*`, `fieller_*` and `mean_omega_*` column; primary columns were unaffected |

None of these may be cited. Their primary estimates from `-a` onward agree with
the current runs to bootstrap noise, so nothing scientific turned on the last two
replacements — but the tables are wrong in the columns named, and a sealed run is
never edited.

## Planned extension, gated on Module 09

Module 10's own scan is null, and the power arithmetic says it was never going
to be otherwise (within-case caudate detects d ~ 1.16 at family-wise
correction; brain smoking-DNAm effects are typically d < 0.3). The one
schizophrenia-relevant question these exposures can answer reframes the
estimand: not "is exposure associated with methylation" but **"is the
risk-variant effect on methylation stable to exposure adjustment"** -- the
stability of an already well-estimated coefficient rather than detection of a
new one.

That would remove the leading alternative explanation for Module 09 (smoking is
72% in cases vs 28% in controls; antipsychotics has zero exposed controls). It
is **not implemented**. Its precondition is now met — Module 09 recorded
accepted runs `scz-AA-{caudate,dlpfc,hippocampus}-20260918` on 2026-09-19 — but
commissioning it remains a PI decision. The design, the mandatory positive control, and the collider
caution are in
`writing-notes/exposure_confounding_of_scz_meqtl_strategy.md`.

## Acceptance gate

`04_apply_gates.R` emits `PASS_EXPLORATORY_COVERAGE` or
`FAIL_EXPLORATORY_COVERAGE`. This is deliberately **not** a scientific success
criterion — a null axis result is a legitimate outcome and must not block
sealing. What it checks is whether the run had the coverage to have said
anything: at least one eligible exposure, at least 500 tested VMRs, at least one
axis model fitted, no forbidden column in any input, and a BH family count that
matches what stage 1 prespecified.

## Accepted runs

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| env-AA-caudate-20260920-a | AA | caudate | vmrset-AA-caudate-937a41979978 | 2026-09-20 | Kynon J.M. Benjamin | PASS_EXPLORATORY_COVERAGE | n=153; 11,251 VMRs; 5/5 coverage criteria; 22/22 chromosomes, 0 excluded/QC-failed/failed/unaccounted. Stage A: 0 FDR-significant VMR x exposure pairs. Stage B: 1 of 9 families survives FDR -- nicotine@all, beta -0.462, q 0.0085, cell_composition_r2 arm -0.440 (p 0.0040). **No percentage here is formally identifiable**: mean_omega_z max 0.78, absolute_p > 0.05 in all 9 families, Fieller unbounded in all 9. **Donor bootstrap inflates the ratio denominator 5.0x-5.8x**, so the bootstrap half of the variance is likely optimistic and these p-values may be too small. Exploratory supplement only (`main_text_retention = NEVER_SUPPLEMENT_ONLY`); the variance-budget limitation (AGENTS.md 7.10) is permanent and a negative gradient may never be read as exposure effects concentrating at weakly controlled VMRs. Caudate is batch-confounded (AGENTS.md 8.1). |
| env-AA-dlpfc-20260920-a | AA | dlpfc | vmrset-AA-dlpfc-856067dfe289 | 2026-09-20 | Kynon J.M. Benjamin | PASS_EXPLORATORY_COVERAGE | n=118; 9,251 VMRs; 5/5 coverage criteria; 22/22 chromosomes, 0 excluded/QC-failed/failed/unaccounted. Stage A: 0 FDR-significant pairs. Stage B: 0 of 5 families survives FDR; nearest is nicotine@schizophrenia (q 0.054). **No percentage is formally identifiable**: mean_omega_z max 0.76, absolute_p > 0.05 in all 5, Fieller unbounded in all 5. **Bootstrap denominator inflation 3.2x-38.5x**, the largest in the module (marital_status@all, the family the retired `relative_scale_min_mean_z` guard was written for, which now simply reports null). Exploratory supplement only; variance-budget limitation applies. |
| env-AA-hippocampus-20260920-a | AA | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-09-20 | Kynon J.M. Benjamin | PASS_EXPLORATORY_COVERAGE | n=117; 9,166 VMRs; 5/5 coverage criteria; 22/22 chromosomes, 0 excluded/QC-failed/failed/unaccounted. Stage A: 12 FDR-significant VMR x exposure pairs (11 nicotine, 1 education) out of ~9,200 tests in 5 families, in the region with the smallest exposed case count -- a supplemental observation, not a finding. Stage B: 2 of 5 families survive FDR -- nicotine@schizophrenia beta -0.491 q 1.1e-5 (arm -0.483, p 1.9e-6) and nicotine@all beta -0.376 q 0.0013 (arm -0.370, p 1.4e-4). **No percentage is formally identifiable**: mean_omega_z max 1.60, never reaching 1.96; absolute_p > 0.05 in all 5; Fieller unbounded in all 5. **Bootstrap denominator inflation 2.0x-4.1x.** Exploratory supplement only; variance-budget limitation applies. |

Accepted by the PI on 2026-09-20. The gate is a **coverage** gate, so what the
acceptance records is that these runs had the coverage to have said something --
not that they found one. Two caveats are written into every row because they
qualify every number in the stage B table and neither is resolved:

1. **No percentage is formally identifiable.** `mean_omega_z` never reaches 1.96
   (maximum 1.60, in hippocampus), so `absolute_p > 0.05` and the Fieller
   interval is unbounded in **every** family in all three regions. The ratio is
   what replicates; the level it is a ratio of is not resolved by these data.
   Never write "49% less exposure-explained variance" without this.
2. **The donor bootstrap inflates the ratio's denominator 2.0×-38.5×.**
   Resampling donors with replacement duplicates donors, which inflates an
   exposure-explained sum of squares, so the bootstrap half of the variance is
   likely optimistic and these p-values may be too small. Settling it needs a
   null plus positive-control simulation that has not been run. The same
   bootstrap construction is in `09b_aging_application`'s accepted runs.

The **variance-budget limitation** (see "Limitation for the discussion") is
permanent and is a Discussion item, not a task: at fixed total variance a higher
genetic share leaves less non-genetic variance for any exposure to move, so a
negative gradient is close to arithmetically forced wherever a real effect
exists. A negative β is never evidence that exposure effects concentrate in
weakly controlled VMRs (AGENTS.md §2.3).

This acceptance unblocks `06_collate_regions.R` without
`--allow-unaccepted-runs`, so `_m/combined/` now carries `citable = TRUE`, and it
unblocks the supplemental `figureS_environmental_axis` in
`11_integrated_manuscript_outputs`. It does **not** move the module out of the
supplement: `main_text_retention = NEVER_SUPPLEMENT_ONLY` stands.

## Contract

`_h/` holds code, `_m/runs/{RUN_ID}/` holds immutable generated output, and
`tests/` is gitignored hand-run smoke checks. Configuration is shared at
`config/environmental.yml`, not duplicated here. See AGENTS.md §5.2 and §7.10.
