# 03_local_snp_prediction — secondary translational endpoint

Answers whether a VMR can be imputed into a new cohort. This is **not** the
primary biological endpoint—that is the continuous local SNP contribution
score from `02_local_genetic_variance`.

**Status: not implemented.** Gated on `02_local_genetic_variance` acceptance (AGENTS.md §6: "No downstream
production run may consume an upstream result until the upstream README records
a passing acceptance gate and immutable run ID").

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

Settings live in `config/prediction.yml`, currently `pi_locked: false`.

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
