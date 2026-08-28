# 03_local_snp_prediction — secondary translational endpoint

Answers whether a VMR can be imputed into a new cohort. This is **not** the
primary biological endpoint—that is the continuous local SNP contribution
score from `02_local_genetic_variance`.

**Status: implemented, smoke-verified, not yet run in production.**
All five stages exist in `_h/` and have been exercised end to end through SLURM
on a 10-VMR smoke run, `lsp-smoke-sched-AA-caudate-20260823` (AA x caudate),
which passed all eight criteria in `04_check_prediction.R` and sealed. That run
establishes only that the code executes, reconciles (10 expected / 9 completed /
1 `qc_failed` / 0 unaccounted / 0 failed) and does not leak (median
`r2_pred_oof` = 0.023, `n_negative` = 5, no clipping). It establishes **nothing
scientific**: 10 loci in one cell is not an estimate of anything, and the run is
marked `smoke_run = TRUE` and is not eligible for the accepted-runs table.

`config/prediction.yml` is `pi_locked: true`: 5 outer x 5 inner folds, 5
repeats, alpha grid {0.1, 0.5, 0.9, 1.0}, `lambda.min`, and a fold-internal
Haseman-Elston screen with a within-fold permutation-calibrated p-value
(1000 permutations, alpha = 0.05).

The `full_data_screen_arm` literature-comparable sensitivity is **deferred, not
overlooked**: `sensitivity.full_data_screen_arm: false` in the locked config.

Production sizing: the smoke measured ~21 s/locus at 199 permutations, so
roughly 30 s/locus at the locked 1000 -- about 96 CPU-hours per cell for 11,530
loci.

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
