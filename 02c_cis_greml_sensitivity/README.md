# 02c_cis_greml_sensitivity — does conventional cis-GREML support Module 02's ordering in the real data?

**Status: accepted 2026-10-06** (`cgs-AA-{caudate,dlpfc,hippocampus}-20261006-a`;
see **Results** and **Accepted runs**). Config PI-locked 2026-10-04 as written,
with two reconciliation-only amendments approved by the PI on 2026-10-06 (see
**Estimator outcomes**).

## Question, and how it differs from 02b

- **02b** asks whether GREML works in this regime at all, on **simulated**
  phenotypes. Answer so far: cis-window REML is unbiased on average, but per
  locus it is imprecise at low h2 (RMSE 0.08–0.12 at true h2 0.05–0.1), and the
  95% CI covers only ~75% of the time at h2 = 0.
- **02c** asks whether conventional cis-GREML, fitted to the **real** VMR
  methylation phenotypes, supports the **existence** and **ordering** of local
  genetic control that Module 02's relative score reports.

## What 02c is not (PI, 2026-10-03)

These are enforced by `assert_cgs_interpretation()`, and every emitted row
carries the flags.

- **It does not replace the Module 02 score.** The score remains the primary
  endpoint. 02c reads it read-only from the accepted run and writes nothing back.
- **It provides no reportable per-VMR absolute h2.** 02b measured per-locus
  error as large as the quantity.
- **It does not define "genetically controlled VMRs" from GREML significance.**
  The LRT is not calibrated per VMR at this n. Estimate, SE, convergence and
  LRT are retained per VMR, and nothing is classified on them.

## Locus model — Module 02's, checked rather than assumed

| | 02c |
|---|---|
| window and SNP QC | ±500 kb, MAF ≥ 0.05, missingness ≤ 0.05, ≥ 100 variants, via `00_shared/locus_io.R::load_observed_locus()`, the reader Module 02 Stage 01 uses |
| SNP set | the QC'd SNP count is compared with Module 02's `num_snps` **per VMR**; the gate fails on any mismatch |
| covariates | the same `age + sex + diagnosis` design matrix, as `--qcovar` |
| GRM | **one** cis-window GRM from exactly those SNPs (`00_shared/gcta.R::grm_from_dosage`, which matches `gcta --make-grm`); no LD-stratified components |
| REML, primary | FUSION's formulation, `--reml --reml-no-constrain --reml-lrt 1`, unconstrained so boundary clipping cannot manufacture a pile at 0 or 1 |
| REML, sensitivity | the same model with Fisher scoring (`--reml-alg 1`). 02b saw GCTA's default AI algorithm diverge in 7 of 9 unconstrained fits at n = 118. Each mode is summarized on its own converged set, and one never fills the other's gaps |
| VMR universe | exactly the accepted Module 02 run's per-VMR table |

## Readouts

All readouts are per mode, on Module 02-eligible VMRs whose fit converged. CIs
come from a delete-one-chromosome weighted block jackknife (the construction
used by Modules 06, 08, 09b and 10).

- **Existence:** the mean cis-GREML estimate across VMRs. 02b found the
  unconstrained estimator unbiased on average, so a mean above zero is evidence
  that local genetic control exists in these data. No single VMR is called.
- **Ordering:** Spearman between the cis-GREML estimate and the Module 02
  score. Concordance with Module 02's HE and BSLMM features is reported
  alongside, descriptively; the score is built from those features, so they are
  not independent of it.
- **Decile profile:** the mean estimate per Module 02 score decile.
- **Selection:** the convergence rate per score decile, so the converged-only
  comparison shows what it conditioned on.
- **LRT:** the fraction with p < 0.05, **descriptive only**.

## Estimator outcomes (amended 2026-10-06)

A fit counts as converged only when GCTA wrote an `.hsq` whose h2 SE is finite
and positive. Three outcomes are estimator outcomes. Each is reconciled as
`qc_failed`, counted per score decile, and kept out of the readouts:
- GCTA's own non-convergence and singular-matrix messages;
- a crash or an "X^t V^-1 X not invertible" error **after** GCTA's iteration
  trace shows a variance component running away. The trace must reach at least
  1e6 x the EM starting Vp. A collinear covariate design stops before
  iterating, so it stays a computational failure;
- an `.hsq` with h2 SE 0 or undefined, which means a degenerate information
  matrix and not a REML optimum. In the superseded round, FUSION fits like
  this reached -128.7.

The two amendments change only how GCTA's exit is classified. The estimator,
modes, window, SNP QC, covariates and summaries are as locked. Both are
opt-in keys in `config/cis_greml_sensitivity.yml`, which records their
counts; 02b's config does not set them.

Per-VMR GREML estimates are noisy (02b: SE about 0.08 at n = 118), which
attenuates any rank correlation with the true ordering. A modest Spearman is
expected even when the orderings agree, and must not be read as disagreement
without that ceiling in view.

## Acceptance gate

`PASS_CIS_GREML_SENSITIVITY_QC` certifies completeness, provenance, and that
02c fitted Module 02's loci and SNPs. It is not a success criterion: a weak or
null ordering is a legitimate result, and the Module 02 score does not change
either way. The six criteria:

1. every unit is reconciled, with zero failures;
2. GCTA is the pinned 1.94.1;
3. the Module 01 and 02 upstreams are still the accepted runs;
4. the SNP set matches Module 02 for every fitted VMR;
5. the interpretation flags are on every row;
6. the primary existence and ordering readouts are reported.

## Pipeline

```
cd 02c_cis_greml_sensitivity/_m && mkdir -p logs
../_h/submit_cis_greml.sh AA caudate      # and dlpfc, hippocampus; SMOKE=1 for a smoke run
Rscript ../_h/06_collate.R                # after acceptance -> _m/combined/
```

Stages: `00_new_run.R` → `01_fit_cis_greml.R` (array, 100 VMRs per task) →
`02_summarize.R` → `03_apply_gate.R` → `04_plot.R` → `05_finalize_run.R`.
Smoke checks are in `tests/test_cgs.R` (gitignored).

## Results (accepted 2026-10-06)

**Conventional cis-GREML supports Module 02's ordering in all three regions.**
Primary mode `fusion_unconstrained`; CIs are delete-one-chromosome jackknife
over 22 blocks.

| region | Spearman with score (FUSION) | Fisher | converged / eligible (FUSION) | mean-estimate CI (existence) |
|---|---|---|---|---|
| caudate | 0.908 [0.895, 0.921] | 0.911 | 10,438 / 11,251 | [0.240, 0.284] |
| DLPFC | 0.914 [0.905, 0.922] | 0.909 | 8,906 / 9,251 | [0.264, 0.294] |
| hippocampus | 0.922 [0.914, 0.929] | 0.915 | 8,786 / 9,166 | [0.288, 0.318] |

- **Ordering.** The mean estimate rises monotonically across the ten score
  deciles in every region and both modes, from -0.007 to -0.025 in decile 1 to
  0.890-0.932 in decile 10.
- **Existence.** The CI excludes zero in 3 of 3 regions. The statement is
  "excludes zero", **never the level**, for two reasons:
  - 02b's unbiasedness was measured on simulated Gaussian phenotypes;
  - decile means are conditioned on a score built from the same data.
  Do not write "on average 26-30% of variance is local" (PI, 2026-10-06).
- **Agreement between estimators on shared data, not replication** (PI,
  2026-10-06).
  - GREML sees Module 02's donors, covariates and cis SNPs: the SNP count
    matches `num_snps` on 11,343 / 9,347 / 9,272 of 11,343 / 9,347 / 9,272
    fitted VMRs.
  - It shows that the score's ordering is not an artifact of the
    elastic-net/joint-model machinery.
  - HE and BSLMM agree with GREML at 0.944-0.960, descriptively. They are
    inputs to the score, so they are expected to sit above it.
- **Selection by convergence, and why it does not drive the ordering.**
  FUSION converges on 60% of caudate's top-decile VMRs, against 91% for Fisher
  scoring, and on 84-86% of every region's bottom decile. Fisher scoring gives
  the same Spearman and decile profile; caudate's decile 10 is 0.915 for
  FUSION and 0.932 for Fisher.
- **Per VMR, nothing is reported or classified.**
  - 49-54% of converged fits have LRT p < 0.05. This is descriptive: 02b
    measured null coverage of 0.72-0.83.
  - 78-81% of estimates are positive.

Reconciliation, FUSION / Fisher `qc_failed` by reason:

| reason | caudate | DLPFC | hippocampus |
|---|---|---|---|
| locus QC (Module 02-ineligible loci; every eligible VMR loaded) | 187 / 187 | 225 / 225 | 225 / 225 |
| information matrix not invertible | 464 / 0 | 214 / 0 | 252 / 0 |
| log-likelihood not converged | 311 / 357 | 131 / 320 | 129 / 288 |
| diverged (crash, or X'V-1X after a runaway trace) | 5 / 0 | 1 / 0 | 4 / 0 |
| degenerate fit, h2 SE 0 or undefined | 49 / 9 | 20 / 4 | 15 / 2 |

Every run has 0 computational failures.

## Accepted runs

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| cgs-AA-caudate-20261006-a | AA | caudate | vmrset-AA-caudate-937a41979978 | 2026-10-06 | Kynon J. M. Benjamin | PASS_CIS_GREML_SENSITIVITY_QC | Built at 7b5b09eb8, git_dirty false, config_cis_greml_sensitivity_sha256 4ae6fbb9... (lock + two reconciliation amendments, PI 2026-10-06); upstreams vmrcat-AA-caudate-20260816, lgv-AA-caudate-rescore-20260913; 23,060 units, 0 failed; SNP set equals Module 02's on all 11,343 fitted VMRs. Ordering: Spearman with score 0.908 [0.895, 0.921] (FUSION, n = 10,438), 0.911 (Fisher). Existence: mean-estimate CI [0.240, 0.284] excludes zero -- not a level. FUSION converges on 60% of the top decile (Fisher 91%, same ordering). Agreement between estimators on shared data, not replication; no per-VMR h2, no significance classes. |
| cgs-AA-dlpfc-20261006-a | AA | dlpfc | vmrset-AA-dlpfc-856067dfe289 | 2026-10-06 | Kynon J. M. Benjamin | PASS_CIS_GREML_SENSITIVITY_QC | Built at 7b5b09eb8, same config; upstreams vmrcat-AA-dlpfc-20260816, lgv-AA-dlpfc-rescore-20260913; 19,144 units, 0 failed; SNP set matched on 9,347/9,347. Spearman 0.914 [0.905, 0.922] (n = 8,906), Fisher 0.909. Existence CI [0.264, 0.294], not a level. Agreement between estimators, not replication; no per-VMR h2, no significance classes. |
| cgs-AA-hippocampus-20261006-a | AA | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-10-06 | Kynon J. M. Benjamin | PASS_CIS_GREML_SENSITIVITY_QC | Built at 7b5b09eb8, same config; upstreams vmrcat-AA-hippocampus-20260816, lgv-AA-hippocampus-rescore-20260913; 18,994 units, 0 failed; SNP set matched on 9,272/9,272. Spearman 0.922 [0.914, 0.929] (n = 8,786), Fisher 0.915. Existence CI [0.288, 0.318], not a level. Agreement between estimators, not replication; no per-VMR h2, no significance classes. |

Accepted by the PI on 2026-10-06 (signed `writing-notes/DRAFT_02c_acceptance_20261006.md`).

## Superseded runs

| run_id | reason |
|---|---|
| cgs-AA-{caudate,dlpfc,hippocampus}-20261004 | Built at ad1326f11 and never sealed. Six FUSION fits (3 / 1 / 2) ended in GCTA exit 1, "X^t V^-1 X not invertible", after AI-REML ran a variance component to 2.9e9-4.2e13 x its EM starting Vp. Fisher converged on all six VMRs. Fixed by amendment 1 (b9f0bb055). |
| cgs-AA-{caudate,dlpfc,hippocampus}-20261006 | Built at b9f0bb055; sealed and passed 6/6. GCTA fits with h2 SE 0 or undefined were counted as converged: 49 / 20 / 15 FUSION and 9 / 4 / 2 Fisher, with FUSION estimates from -128.7 to +5.2. One of them moved hippocampus's existence mean from 0.303 to 0.288 and made a spurious decile-3 dip. Fixed by amendment 2 (7dd7fbc88). In the -a round, completed fell by exactly 58 / 24 / 17, the degenerate counts, and nothing else changed. |
| cgs-AA-dlpfc-smoke-20261003 | Smoke run; never intended for acceptance. |

**Process note.** At about 08:46 on 2026-10-06, while the -a fit arrays were
running, the main working tree was checked out to another branch for under a
minute. None of the runs' inputs differed:
- `00_shared/` is identical on both branches;
- fits read phenotypes from the Module 01 run, not from `config/paths.yml`;
- the one job that started in that window reads only run-directory files.

`git_commit` and `git_dirty` were recorded at submission. Branch work is done
in worktrees from now on.

## Contract

This module follows AGENTS.md §5.2. Its configuration is
`config/cis_greml_sensitivity.yml`. It depends on 01 and 02 and reads both
read-only.
