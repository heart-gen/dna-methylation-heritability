# 10_environmental_exploratory — exploratory environmental-factor analysis

**Status: implemented, no run submitted. Supplement only, permanently.**

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
| **tobacco** (smoking \| nicotine) | 76 / 0.000 | 58 / 0.000 | 57 / 0.000 | all three |
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

#### The stratified floor is 15, and was chosen after seeing the counts

PI decision 2026-09-10. The pooled floor stays at 20 and was locked before any
count was seen. The stratified floor is 15, and it was lowered from 20 **after**
the within-schizophrenia counts were computed — it admits `antipsychotics` (18
exposed of 65) and `tobacco` (19), both of which missed 20 narrowly.

A stratum halves the sample by construction, so a lower floor is defensible on
its own terms. But the sequence is what it is, and `config/environmental.yml`
records `decided_after_seeing_counts: true` so a reader can weigh it. It does
not admit `antipsychotics` in dlpfc or hippocampus (11), `education` (12) or
`marital_status` (7) anywhere.

Eligible exposure × stratum pairs, AA caudate:

| stratum | eligible |
|---|---|
| `all` | tobacco, any_substance_nontobacco, education, marital_status |
| `schizophrenia` | tobacco (19), antipsychotics (18), any_trauma_hx (22) |

In dlpfc and hippocampus the schizophrenia stratum is n=48 and nothing clears
the gate; those regions run the pooled stratum only.

### Drug-use composites

Regrouping the substance columns into pharmacologically coherent classes does
not recover power: opioids (codeine|morphine|fentanyl) reach 7/6/6, and
psychostimulants (cocaine|amphetamines) 9/6/6, across the three AA regions. Both
are recorded in `exposures.retired` with their counts so the attempt is not
silently repeated. Aggregation **across** mechanism is what crosses the gate, and
that yields two changes:

- **`tobacco`** replaces `smoking` and `nicotine` as separate exposures. They are
  the same exposure recorded two ways and overlap heavily; testing both reports
  one signal twice. They remain as descriptive rows outside the BH families.
- **`any_substance_nontobacco`** is a new eligible exposure in AA caudate only.
  It is a **burden indicator, not a pharmacological class** — an opioid-positive
  and an alcohol-positive donor share no mechanism — and no statement about any
  single drug is licensed by a result on it. dlpfc and hippocampus return 19
  against a gate of 20; `exposure-eligibility.tsv` flags that as `near_gate`
  rather than presenting a one-donor miss as a clean exclusion.

An open item before the first production run: confirm against the LIBD donor
metadata whether each source column is toxicology at death or lifetime history.
`tobacco` currently unions one of each. If the windows differ,
`any_substance_nontobacco` stays toxicology-only and the `tobacco` union is
revisited.

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

    cd 10_environmental_exploratory/_m && mkdir -p logs
    ../_h/submit_environmental.sh AA caudate

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

1. **Primary, threshold-free:** `-log10(p) ~ score_z + vmr_length + cpg_count +
   cpg_density + gc_content + mappability + mean_methylation +
   methylation_variance`. Continuous on both sides, so the headline statement
   does not depend on where the FDR cut lands.
2. Wilcoxon of `score_z`, FDR-significant vs not.
3. Logistic on the same indicator with the same covariates.

(2) and (3) are secondary — they need a threshold — and are skipped with a
recorded reason when either group holds fewer than 10 VMRs. Technical covariates
are joined from Module 04's `vmr-features.tsv` rather than recomputed, so GC
content and mappability are identical to the values Module 04 published.

## Interpretation constraints

- A null or negative axis result is **not** evidence that a VMR is
  environmentally determined (AGENTS.md §2.3).
- No heritable/non-heritable groups, no absolute PVE, no legacy metric.
  `forbidden_columns` is checked against every upstream table at runtime.
- **Region-specific only.** Caudate is sequencing batch 3 and region is
  perfectly confounded with batch (AGENTS.md §8.1), so no cross-region exposure
  contrast is emitted and none is licensed.
- `antipsychotics` does not clear the gate in any region, so it is not tested
  here. It stays declared because it is downstream of diagnosis, which is in the
  locus model — a collider — and that flag must survive if the missingness
  ceiling is ever revisited.
- `any_substance_nontobacco` is a burden indicator across unrelated mechanisms.
- Supplement only, always. The decision row carries
  `main_text_retention = NEVER_SUPPLEMENT_ONLY`.

## First end-to-end run (smoke, 2026-09-10)

`env-smoke2-AA-caudate-20260910` ran all stages on the accepted upstreams. It is
a **smoke run and must not be cited**: it was driven stage-by-stage rather than
through the SLURM chain, so it does not satisfy AGENTS.md §9. It is recorded
because it establishes that the pipeline works, and because its shape is the
most informative thing the module has produced.

22/22 chromosomes reconciled; 11,341 VMRs per family; 7 exposure × stratum
families; `PASS_EXPLORATORY_COVERAGE`.

**Stage A is null in both strata.** Zero FDR-significant VMRs in all seven
families. The smallest p anywhere was 1.3 × 10⁻⁵ (tobacco, pooled) across 11,341
tests — nothing close to surviving BH.

**Stage B splits by stratum, and that is the finding:**

| exposure | stratum | beta | p |
|---|---|---|---|
| tobacco | all | −0.0328 | 3.0 × 10⁻⁸ |
| marital_status | all | 0.0181 | 1.9 × 10⁻⁴ |
| any_substance_nontobacco | all | 0.0164 | 8.1 × 10⁻⁴ |
| education | all | 0.0088 | 0.079 |
| tobacco | schizophrenia | 0.0021 | **0.67** |
| antipsychotics | schizophrenia | 0.0052 | **0.29** |
| any_trauma_hx | schizophrenia | 0.0039 | **0.43** |

The pooled gradient is nominally strong and **vanishes inside the case
stratum** — tobacco goes from 3 × 10⁻⁸ to 0.67, with the coefficient collapsing
to a fifteenth of its pooled size.

**Read this as a warning, not a result.** Three things have to be held together:

1. Stage A found nothing, so Stage B is describing the shape of a p-value
   distribution that is entirely noise. A slope in null p-values is not an
   exposure effect.
2. The grouped tests could not run at all (`fewer_than_10_vmrs_in_a_group`), so
   there is no second line of evidence for any row above.
3. The stratum contrast has a mundane competing explanation — n=65 against
   n=153 — that this design cannot separate from the interesting one, that the
   pooled gradient tracks something correlated with diagnosis rather than with
   exposure.

Taken together, the pooled `p = 3 × 10⁻⁸` is far more likely to be a residual
technical gradient along the score than evidence that exposure effects
concentrate in low-genetic-control VMRs. That inference is exactly what
AGENTS.md §2.3 forbids, and it is why the acceptance gate is a coverage check
rather than a success criterion.

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

None. No production run has been submitted.

## Contract

`_h/` holds code, `_m/runs/{RUN_ID}/` holds immutable generated output, and
`tests/` is gitignored hand-run smoke checks. Configuration is shared at
`config/environmental.yml`, not duplicated here. See AGENTS.md §5.2 and §7.10.
