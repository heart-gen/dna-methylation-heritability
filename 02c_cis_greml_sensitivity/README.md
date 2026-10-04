# 02c_cis_greml_sensitivity — does conventional cis-GREML support Module 02's ordering in the real data?

**Status: implemented 2026-10-03; config PI-locked 2026-10-04 as written; no accepted run yet.**

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

## Accepted runs

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
|---|---|---|---|---|---|---|---|
| _(none)_ | | | | | | | |

## Contract

This module follows AGENTS.md §5.2. Its configuration is
`config/cis_greml_sensitivity.yml`. It depends on 01 and 02 and reads both
read-only.
