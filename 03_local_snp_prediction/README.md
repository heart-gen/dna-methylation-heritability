# 03_local_snp_prediction — secondary translational endpoint

Answers whether a VMR can be imputed into a new cohort. This is **not** the
primary biological endpoint—that is the continuous local SNP contribution
score from `02_local_genetic_variance`.

**Status: run in production; all three AA cells accepted 2026-08-28.**
See the **Accepted runs** table below for `lsp-AA-{caudate,dlpfc,hippocampus}-20260825`.
Each reconciled with zero computational failures and zero unaccounted loci, and
each carries the same headline result: median `r2_pred_oof` ~ 0 (-5.7e-07,
-2.1e-06, -5.2e-07), mean 0.17-0.20, 44-47% of loci positive. **These values
support relative ranking of loci only**, not a claim that any individual VMR is
imputable into a new cohort. Post-hoc QC on the predictability distribution
lives in `_h/06_qc_predictability_artifact.R` (see below).

The code was first exercised end to end through SLURM on a 10-VMR smoke run,
`lsp-smoke-sched-AA-caudate-20260823` (AA x caudate), which passed all eight
criteria in `04_check_prediction.R` and sealed. That run established only that
the code executes and does not leak; it is marked `smoke_run = TRUE` and is not
eligible for the accepted-runs table.

### Known caveat: `calibration_slope` in the accepted runs

`_h/03_combine_oof.R` guards the per-locus calibration slope against
numerically degenerate loci -- a locus screened in for a handful of folds and
given the prespecified null prediction elsewhere has a predictor spread of order
1e-16, and dividing by it yields a slope of order 1e16. That guard landed in
`ed6abe0b8` (2026-08-28), **after** the three accepted runs were produced, so
their per-locus tables have no `calibration_slope_defined` column and contain
values up to 5.25e17 (135 such loci in caudate, 110 in dlpfc, 89 in
hippocampus, taking |slope| > 1e6).

This does not affect any reported quantity. The run-level statistic is
`median_calibration_slope`, which is robust to the degenerate tail and reads
1.087 / 1.069 / 1.107 across the three regions. No downstream module reads the
per-locus column. **Do not take a mean over `calibration_slope` in these
runs**; regenerate under the current code first, or restrict to loci that
would clear the tolerance.

`config/prediction.yml` is `pi_locked: true`: 5 outer x 5 inner folds, 5
repeats, alpha grid {0.1, 0.5, 0.9, 1.0}, `lambda.min`, and a fold-internal
Haseman-Elston screen with a within-fold permutation-calibrated p-value
(1000 permutations, alpha = 0.05).

The `full_data_screen_arm` literature-comparable sensitivity is **deferred, not
overlooked**: `sensitivity.full_data_screen_arm: false` in the locked config.

Production sizing: the smoke measured ~21 s/locus at 199 permutations, so
roughly 30 s/locus at the locked 1000 -- about 96 CPU-hours per cell for 11,530
loci.

## Post-hoc QC: is the predictability distribution an artifact?

`_h/06_qc_predictability_artifact.R`, launched by `_h/step_6_qc_predictability.sh`,
runs three independent checks on an accepted run and writes to `_m/qc/{run_id}/`
(gitignored -- regenerable from the sealed run plus `_h/`). It deliberately runs
from the live `_h/` rather than a run snapshot, because it is post-hoc analysis
*of* a sealed run, not part of it.

The motivating observation was a right tail of loci with `r2_pred_oof` > 0.9,
markedly heavier in caudate. The checks:

- **A. Annotation.** Length-matched overlap comparisons plus logistic regression
  adjusted for `log_len` and `median_n_variants`. The segmental-duplication /
  CNV hypothesis was **refuted**: high-r2 loci show *less* segdup overlap than
  length-matched controls (0.043 vs 0.173, p = 3.8e-08), zero blacklist overlap,
  and no mappability deficit. What does hold is a monotone LINE/L1 gradient
  across r2 bins -- caudate 0.048 / 0.127 / 0.414 / 0.706, reproducing in dlpfc
  and hippocampus -- with adjusted OR 3.29 (p = 5e-17). **This is the finding
  that makes `line_l1_frac` x `r2_pred_oof` near-circular in Module 04**, and is
  why that cell is descriptive-only there.
- **B. Relatedness.** plink2 KING over 175,664 LD-pruned autosomal variants,
  11,628 pairs across 153 donors: **maximum kinship 0.0425**, below the
  third-degree threshold (0.0442). No cryptic relatedness leaks across folds.
- **C. Cross-region concordance and a threshold ladder.** Caudate's 347 high-r2
  loci reach median r2 0.829 in DLPFC with 99.2% above 0.5; Spearman 0.784 over
  4,047 shared loci. The caudate excess shrinks smoothly with the cutoff (25x at
  0.90, 5.4x at 0.85, 2.7x at 0.80, 1.5x at 0.70, 1.1x at 0.50), which is a
  uniformly shifted distribution crossing an arbitrary threshold -- consistent
  with caudate's larger donor count (153 vs 118) -- not a discrete artifact
  population.

Conclusion: the tail is not an artifact of repeats-driven mismapping,
relatedness, or region-specific pathology. It is not, however, evidence that
those loci are individually imputable; the relative-ranking restriction above
still applies.

## Migrating from

`local-snp-prediction/`, plus the untracked `local-snp-prediction/oof_diagnostic/` that established defect E1.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## Why the legacy module cannot be reused

`local-snp-prediction/*/01.elastic-net.R` reports

```r
r_squared_cv <- cor(pheno_scaled, predict(final_model, G_clumped[], s = "lambda.min"))^2
```

This is an **in-sample** fit: `cv.glmnet` used folds only to pick lambda, and SNP
clumping was supervised on all donors. Reported values run ~0.85; the honest
held-out value is ~0.01, and `r_squared_cv` correlates *negatively* with `r2_oof`
(rho ~ -0.05 to -0.07). It is invalid as an input to any v2 model, table, or
figure (AGENTS.md §3).

## Required design

The primary estimate must be **end-to-end out-of-fold**. Within each outer split,
held-out donors must not influence phenotype residualization/centering/scaling,
genotype imputation/scaling/MAF/missingness/zero-variance filters,
phenotype-informed SNP selection, locus-level cis screening, or any tuning.

A VMR that fails its fold-internal screen receives the prespecified null
prediction for held-out donors — failed folds are not dropped.

Required metrics, all from pooled OOF predictions: `r2_pred_oof = 1 - SSE/SST`,
`cor2_oof`, RMSE, MAE, calibration intercept and slope, screening pass frequency,
predictions per donor, and fold/repeat diagnostics. **Negative `r2_pred_oof` is
retained**; `cor2_oof` is never substituted when `r2_pred_oof` is unfavorable.

Do not use the current ordinary-OLS Haseman-Elston p-value as the production
screen without donor-robust or permutation calibration. Lock the screen and
multiple-testing rule before viewing final results.

A full-data-screen-then-nested-OOF arm may be included as a literature-comparable
sensitivity, explicitly labeled.

Settings live in `config/prediction.yml`, now `pi_locked: true` (see Status above).

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.

## Accepted runs

| run_id                      | cohort | region      | vmr_set_id                         | accepted_on | accepted_by         | decision               | notes                                                                                                                               |
|-----------------------------|--------|-------------|------------------------------------|-------------|---------------------|------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
| lsp-AA-caudate-20260825     | AA     | caudate     | vmrset-AA-caudate-937a41979978     | 2026-08-28  | Kynon J.M. Benjamin | PASS_OOF_PREDICTION_QC | 11530 expected / 11343 scored / 187 QC-failed / 0 failed; median r2_pred_oof ~0, mean 0.204, 45.1% positive; relative rank use only |
| lsp-AA-dlpfc-20260825       | AA     | dlpfc       | vmrset-AA-dlpfc-856067dfe289       | 2026-08-28  | Kynon J.M. Benjamin | PASS_OOF_PREDICTION_QC | 9572 expected / 9347 scored / 225 QC-failed / 0 failed; median r2_pred_oof ~0, mean 0.174, 44.4% positive; relative rank use only   |
| lsp-AA-hippocampus-20260825 | AA     | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-08-28  | Kynon J.M. Benjamin | PASS_OOF_PREDICTION_QC | 9497 expected / 9272 scored / 225 QC-failed / 0 failed; median r2_pred_oof ~0, mean 0.188, 46.5% positive; relative rank use only   |
