# 02b_greml_simulation_benchmark — can REML recover absolute local h2 at this cohort's n?

**Status: implemented and PI-locked 2026-10-03; no accepted run yet.**
`config/greml_benchmark.yml` was locked as proposed: the scenario grid below,
unconstrained Fisher scoring primary, v1's constrained AI-REML secondary.

## Question

How well does a classical REML estimator (GCTA-GREML) recover **absolute**
local SNP heritability at this cohort's donor counts and real cis-window LD?

AGENTS.md 7.2 retired absolute locus PVE after the frozen elastic-net joint
model failed 6 of 14 absolute-PVE gates while passing the ordering gate. That
evidence concerns one estimator. This module asks the question of an estimator
that shares nothing with it. It is a port of the legacy `simulation-analysis/gcta/`
benchmark (PI decision 2026-10-03).

The answer is not assumed in either direction. If REML recovers absolute local h2
well on simulated phenotypes, that narrows 7.2's retirement to the estimator it
was measured on; it does not by itself license an absolute PVE for any observed
VMR, which this module never estimates (see **What it does not do**). That
question, if it is to be asked, is a PI decision and a new analysis.

## What it does not do

- **It never estimates an observed methylation locus.** Every phenotype is
  simulated, so no observed-locus PVE can come out of it (AGENTS.md §3, §7.2).
  Every row carries `simulated_phenotypes_only = TRUE` and
  `absolute_pve_interpretation_allowed_for_observed_loci = FALSE`.
- **It does not reopen Module 02.** It reads nothing from and writes nothing to
  `02_local_genetic_variance/`. Module 02's runs, score and
  `PASS_RELATIVE_GENETIC_CONTROL_FAIL_ABSOLUTE_LOCUS_PVE` decision are untouched,
  and nothing here can change a Module 02 rank.
- **It is not the observed-score cross-check.** Correlating a second estimator
  against the accepted `local_snp_contribution_score` (v1
  `simulation-analysis/comparison/correlation/`, D4) was closed by PI decision
  on 2026-10-03, not migrated.
- **No classes.** Recovery is reported continuously. The legacy correlation
  script's Heritable / Non-heritable / Low prediction split, defined on a 0.1
  cut, is not ported (AGENTS.md §2.3, §3).

## Design

Two arms, read on the same REML modes.

### Arm 1 — faithful v1 port, `ld_model = ar1_out_of_regime`

One run, `greml-sim-ar1-{date}`, over every simulated sample size in
`inputs/simulated-data/_m/sim_{100,150,200,250,500,1000,5000,10000}_indiv`
(1,000 phenotypes each). The v1 steps are ported unchanged in substance:

| v1 | v2 |
|---|---|
| `step_1a.sh`: `--ld-score-region 200` per chromosome | `_h/02_arm1_ldscore.R` |
| `step_1b.sh` + `01.stratify_LD.R`: quartile LD groups at `summary()` Q1/median/Q3, four GRMs | `_h/03_arm1_stratify_grm.R` (stratification verified identical in `tests/`) |
| `step_1c.sh`: `--reml --mgrm --mpheno k`, no covariates | `_h/04_arm1_reml.R` |
| `02.summary.py`: melt `.hsq`, report `Sum of V(G)/Vp` | `_h/greml_functions.R::parse_hsq()` |

What changed:

- paths come from config;
- seeds and run IDs come from `00_shared/runid.R`;
- every unit is reconciled;
- outputs go to immutable runs;
- the phenotype file is checked against `pheno_1..pheno_K` and passed without
  its header row. v1 passed the header straight to `--pheno`, where GCTA read
  it as a donor named `FID`.

**Arm 1 is out of regime, for two reasons, and its numbers describe the AR(1)
design only** (`interpretation.arm1_transfers_to_cohort = false`):

1. Genotypes are AR(1) Toeplitz LD
   (`inputs/simulated-data/_h/01.simulated-data.py:31`). AR(1) LD decays
   exponentially, and its effective SNP count runs far above real cis-windows.
   Module 02's own AR(1) grid gave low-PVE conclusions that reversed on real
   genotypes.
2. GREML-LDMS fits four **genome-wide** LD-stratified GRMs to a signal confined
   to a 1-Mb window.

The v1 summary table shows how often v1 REML returned nothing at small n. It
holds 1 row at n = 150 and 1 at n = 200, against 570 at n = 10,000. v1 kept no
record of why. Here every non-converged fit is reconciled as `qc_failed` with
GCTA's message, and the rate is reported per cell.

### Arm 2 — observed regime, `ld_model = observed_AA_cis`

One run per region, `greml-AA-{region}-{date}`. Each run reads real AA
cis-window genotypes from the accepted Module 01 catalog run. Only the
phenotype is simulated, through the same two shared helpers Module 02's
observed-regime grid used:

- `00_shared/locus_io.R::load_observed_locus()`, which gives the locked cis
  window (±500 kb), SNP QC (MAF ≥ 0.05, missingness ≤ 0.05, ≥ 100 variants),
  donor alignment, `design_n` (153 / 118 / 117), and the
  `age + sex + diagnosis (+ snpPCs)` covariates;
- `00_shared/locus_io.R::simulate_phenotype_on_observed_genotype()`.

The scenario grid uses Module 02's locked observed-regime axes, so the two
estimators are characterised on the same scenarios:

- h2 ∈ {0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.95};
- three architectures;
- 2 replicates per locus and cell.

Loci are drawn per region by quartile of pre-QC cis SNP count, 50 per quartile.
The draw uses a fixed seed namespace, so a superseding run draws the same loci.

The GRM is written directly from the dosage matrix `load_observed_locus()`
returns (`greml_functions.R::grm_from_dosage()`). It reproduces
`gcta64 --make-grm` on the same SNPs and donors to float32 precision: max
|diff| 1.0e-7 on a 3,428-SNP DLPFC locus with 1,930 missing calls (repeated in
`tests/`). REML adjusts for the locked covariates through `--qcovar`.

### REML modes (PI-locked 2026-10-03)

| mode | role | GCTA flags | why |
|---|---|---|---|
| `unconstrained_fisher` | **primary** | `--reml-no-constrain --reml-alg 1` | unbiased estimator allowed below 0; keeps low-h2 information, as Module 02's unbounded rank does |
| `constrained_ai` | secondary, the v1 setting | (GCTA default) | what v1 ran; piles estimates onto the 0 bound |

These modes come from a 2026-10-03 prototype on one real DLPFC locus
(n = 118, 3 h2 values × 3 replicates):

- unconstrained **AI**-REML diverged in 7 of 9 fits (V(e) went negative and the
  information matrix became singular), which is why it is not offered;
- unconstrained Fisher scoring converged in 9 of 9;
- constrained AI-REML converged in 9 of 9.

In the smoke tests, the primary mode converged in 19 of 20 null fits, with mean
estimate 0.014.

## Why v1's GREML failed at small n, and why the cis-window arm does not

The v1 benchmark (`simulation-analysis/gcta/`) was the staff analysis that found
GREML failing or performing very badly at limited sample sizes. That result is
real, and arm 1 reproduces it. The v1 summary, joined to the simulation truth:

| n | fits that returned anything (of 1,000) | median SE of h2 | estimates at 0 or 1 | Spearman(truth, estimate) |
|---|---|---|---|---|
| 150 | 1 | 15.0 | 100% | -- |
| 250 | 2 | 9.6 | 50% | -- |
| 500 | 18 | 4.7 | 56% | 0.02 |
| 1,000 | 94 | 2.4 | 54% | -0.20 |
| 5,000 | 515 | 0.46 | 9% | 0.35 |
| 10,000 | 570 | 0.23 | 2% | 0.57 |

The missing fits are constrained REML stopping at the 0/1 bound. The arm 1
smoke run at n = 100 and 150 reproduced this: 0 of 20 constrained fits finished,
and unconstrained estimates ran from -42 to +54 with SEs of 8-29.

**The cause is information, not the estimator.** GREML's information comes from
the spread of the GRM's off-diagonal entries. For unrelated people,
SE(h2) is approximately sqrt(2) / (n x SD of the off-diagonal GRM entries).

- **Genome-wide GRM.** Tens of thousands of effectively independent segments
  make the off-diagonal SD tiny (about 0.0045), so SE is about 316/n: roughly 3
  at n = 100 and 0.3 at n = 1,000. Splitting it into four LD-stratified
  components makes this worse.
- **Cis-window GRM.** For one real DLPFC locus (3,428 SNPs, n = 118), LD leaves
  about 48 effective segments. The off-diagonal SD is 0.145, so the predicted SE
  is about 0.08. Arm 2's measured RMSE (0.08-0.12 at true h2 0.05-0.1) matches.

v1 had two further handicaps:

- It was misspecified. The GRM spreads the genetic signal over 7.7M SNPs while
  the signal comes from 1-5 SNPs in one 1-Mb window.
- Its AR(1) LD is out of regime for real cis windows.

**What arm 2 does and does not show.**

- The cis-window estimator is unbiased on average, and it orders a wide range of
  true h2 well.
- Per locus, though, at the h2 most VMRs sit at, its error is as large as the
  quantity, and its CI under-covers at h2 = 0.
- Arm 2 is also a best case:
  - the GRM is built from exactly the causal SNPs' window;
  - the noise is Gaussian;
  - the realized h2 is fixed exactly;
  - there is no measurement error;
  - covariates are exactly linear.

So v1's conclusion still holds for per-VMR absolute values. It does not extend
to cis-window REML as such, which is how FUSION-style analyses screen
cis-heritability at a few hundred samples. Whether conventional cis-GREML
supports Module 02's ordering on the **real** phenotypes is the separate question
asked by `02c_cis_greml_sensitivity`.

## Metrics

All metrics are against the **realized** simulated h2, per cell (arm 1: n × mode;
arm 2: mode × architecture × h2):

- bias, RMSE, 95% CI coverage;
- boundary rate (constrained) or rate outside [0, 1] (unconstrained);
- Spearman with a cluster-bootstrap CI (clusters are loci in arm 2 and
  phenotypes in arm 1);
- the REML estimation-failure rate.

## Acceptance gate

`PASS_GREML_BENCHMARK_QC` certifies **completeness and provenance, not good
recovery.** Poor recovery at small n is a legitimate result — it is the result
this module exists to measure — and does not block sealing. The six criteria,
in `summary/gate-criteria.tsv`:

1. every unit is reconciled, with zero failed, unaccounted or unexpected units;
2. GCTA is the pinned 1.94.1, both when the run opened and at the gate;
3. arm 1 inputs are checksummed, or the arm 2 Module 01 upstream is still
   accepted;
4. primary-mode fits exist at both ends of the truth range;
5. the simulated-only flags are on every row;
6. primary-mode metrics and the failure rate are reported in every cell.

## Interpretation constraints

- Arm 2 describes REML **on simulated phenotypes** at this cohort's n and LD. It
  says how well absolute local h2 *could* be recovered. It is never an estimate
  for any VMR.
- Do not write a percentage of methylation variance explained for any locus or
  class from this module.
- Arm 1 does not transfer to the cohort. Quote it only as the reproduction of
  the v1 benchmark, with both out-of-regime reasons.
- Nothing here reorders, re-scores or validates a Module 02 locus.

## Pipeline

```
cd 02b_greml_simulation_benchmark/_m && mkdir -p logs
../_h/submit_greml_benchmark.sh ar1
../_h/submit_greml_benchmark.sh observed AA caudate      # and dlpfc, hippocampus
# smoke: SMOKE=1 ...   plan only: DRY_RUN=1 ...
# after acceptance:
Rscript ../_h/10_collate.R                                # -> _m/combined/
```

| stage | script | unit |
|---|---|---|
| open | `00_new_run.R` | run, task manifests, locus draw |
| arm 1 | `01_arm1_inputs.R`, `02_arm1_ldscore.R`, `03_arm1_stratify_grm.R`, `04_arm1_reml.R` | inputs; n × chromosome; n; phenotype chunk |
| arm 2 | `05_arm2_greml.R` | locus chunk |
| both | `06_summarize.R` → `07_apply_gate.R` → `08_plot.R` → `09_finalize_run.R` | run |
| collate | `10_collate.R` | accepted runs |

Environment:

- GCTA 1.94.1 from bioconda, in `/gpfs/projects/p32505/opt/envs/genomics`.
  It is called by absolute path from the `epigenomics` R env, and no
  `module load` is used.
- Smoke tests are in `tests/test_greml_functions.R` (gitignored, §5.2).

## Accepted runs

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| _(none)_ | | | | | | | |

## Migrating from

- `simulation-analysis/gcta/` → arm 1, with v1 parameters, plus arm 2.
- `simulation-analysis/comparison/correlation/` → **closed**. Its
  simulated-truth GREML correlation is subsumed by this module's continuous
  metrics. Its `h2_category` split is withdrawn, and the observed-score
  cross-check is not migrated (PI 2026-10-03).

See `MIGRATION_MANIFEST.tsv`.

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, and `tests/` holds
gitignored smoke checks. Its shared configuration is `config/greml_benchmark.yml`.
