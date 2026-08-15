# Phase 2: VMR-level aggregation of CpG meQTL evidence

## Question

Does CpG meQTL burden increase with continuous VMR local genetic predictability?

## Run

```bash
cd meqtl-validation/02_vmr_meqtl_burden/_m
mkdir -p logs
sbatch ../_h/step_1.sh   # aggregate (qval≤0.05 significance; joins coverage + technical annotations)
# Then run Phase 6 step_2_tech_joins.sh in join-only mode.
sbatch ../../07_repeat_mappability_sensitivity/_h/step_2_tech_joins.sh
sbatch ../_h/step_2.sh   # overdispersion-robust models + matched-pair inference
```

Primary outcome: proportion of CpGs with `qval ≤ 0.05`.  
Primary predictor: continuous VMR local genetic predictability (`h2_unscaled`).

**Phase 1 lock:** primary calls are **M3a** (`../01_cpg_meqtl_mapping/PHASE1_LOCK_DECISION.md`).  
Prior M0 Phase 2 outputs archived at `_m/M0_archive_20260730/`.

## Repair-v2 status

The existing `_m/` tables below are historical pre-repair outputs and must not be
used for claims. Schema-v2 aggregation counts only CpGs present in the TensorQTL
tested-results table, preserves untested VMRs for audit, and fits overdispersion-
robust models plus propensity-score matched-pair randomization tests. Downstream
modules reject burden tables with `analysis_schema_version < 2`.

## Calibrated-estimator comparison workflow

The accepted simulation-calibrated local SNP variance estimates are evaluated as
an additive sensitivity arm; they do not overwrite the manuscript's legacy
local-predictability tables. The locked comparison uses `observed-AA-v4`, exact
one-to-one VMR coordinates, M3a CpG-level meQTL calls, the repair-v2 denominator,
and the same complete VMR set for calibrated and legacy regressions.

Submit the full immutable workflow from the repository root:

```bash
bash meqtl-validation/02_vmr_meqtl_burden/_h/submit_predictor_comparison_workflow.sh \
  calibrated-vs-legacy-AA-v1
```

The submission helper explicitly uses `/projects/p32505/opt/envs/genomics`,
verifies required Python packages and the accepted calibrated-run QC, snapshots
code/configuration and accepted calibrated inputs, records input SHA-256 values,
and submits this dependency chain:

```text
repair-v2 aggregation (3-region array)
  -> technical-coordinate join
  -> identical-set legacy/calibrated models (3-region array)
  -> combined FDR, promotion gates, and QC figure
```

Run-scoped results are written under:

`_m/predictor-comparison-runs/<RUN_ID>/`

Important outputs are:

- `regions/<region>/predictor_bridge.tsv.gz`: exact joined analysis table;
- `combined/burden_model_comparison.tsv`: standardized head-to-head effects;
- `combined/matched_comparison.tsv`: matched positive/nonpositive contrasts;
- `combined/predictor_concordance.tsv`: estimator concordance diagnostics;
- `combined/comparison_acceptance.tsv`: prespecified Phase 2 gates;
- `combined/promotion_status.tsv`: explicit manuscript-promotion status;
- `figures/predictor-comparison.{pdf,svg,png}`: manuscript-facing QC figure;
- `provenance/`: environment, commit, code/config, inputs, checksums, and jobs.

Interpretation is hierarchical. The calibrated estimate may be described as
local SNP-explained methylation variance only if the calibrated Phase 2 gate
passes. Full replacement of the manuscript axis additionally requires the
prespecified LINE/L1 and repressive-chromatin sensitivity analyses, which are
deliberately marked pending by this workflow. Until then, report the calibrated
analysis beside the legacy aggregate local SNP contribution.

## Calibrated-estimator promotion follow-up

After `calibrated-vs-legacy-AA-v1` passes its Phase 2 burden gate, run the
prespecified promotion addendum:

```bash
bash meqtl-validation/02_vmr_meqtl_burden/_h/submit_predictor_followup_workflow.sh \
  calibrated-promotion-AA-v1
```

This immutable workflow does not modify the accepted comparison run. It:

1. reaggregates only locked M3a-tested CpGs into the pre-existing
   all-individual VMR coordinates (untested CpGs are never negatives);
2. repeats burden models after excluding calibrated upper-boundary hits;
3. uses propensity overlap weighting as the prespecified alternative to the
   narrowly imbalanced hippocampal matched analysis;
4. compares calibrated and legacy predictors for LINE/L1, H3K9me3, and
   quiescent-chromatin enrichment under high-mappability, SNP-proximity, and
   segmental-duplication restrictions;
5. tests independent Jaffe DLPFC meQTL support while retaining Schulz
   hippocampal overlap as positive-only descriptive evidence;
6. reuses existing expression, PSI, and ABC association tables without
   rerunning transcriptome scans;
7. evaluates locked promotion gates and exports a manuscript-ready figure.

Outputs are written under
`_m/predictor-followup-runs/<RUN_ID>/{coordinate_aligned,sensitivities,annotation,orthogonal,combined,figures,provenance}`.
The exact-intersection AA-VMR analysis remains primary because locked M3a did
not test every CpG in the all-individual VMR universe. Full estimator promotion
requires every criterion in `combined/followup_acceptance.tsv` to pass.

### Completed locked run (2026-08-08)

`calibrated-promotion-AA-v1` completed successfully under SLURM. The calibrated
estimate passed the internal Phase 2, coordinate-aligned, upper-boundary,
overlap-weighting, and repeat/repressive-chromatin gates. Existing expression
associations were also positive in all three regions. The Jaffe DLPFC join
resolved 1,050 common assayed VMRs, but all 1,050 were supported; the outcome
therefore has no variation and the external gradient is not estimable. The
positive-only Schulz resource similarly cannot provide negative controls.

Under the prespecified gate, full replacement is **not accepted**. Retain
`h2_unscaled` as the locked manuscript axis and report `h2_en_calibrated` as a
strong, directionally consistent sensitivity result. Do not recode non-assayed
VMRs as external negatives or weaken the external gate after seeing the result.
The audited decision is in
`_m/predictor-followup-runs/calibrated-promotion-AA-v1/combined/promotion_status.tsv`.
The modular Manubot-ready Methods, Results, figure/table, reproducibility, and
limitations summary is `ANALYSIS_SUMMARY_CALIBRATED_PROMOTION.md`.

## Historical results (invalidated pending rerun; 2026-08-01 M3a calls)

Binomial GLM: meQTL burden ~ z(local_predictability) [+ covariates]. Positive coef = higher predictability → higher meQTL burden.

| Region | N VMRs (with pred.) | Unadj. coef | Tech-adj. coef* | Matched Δ (hi−lo) | Matched perm. p |
|---|---:|---:|---:|---:|---:|
| Caudate | 9479 | 3.35 | 3.20 | +0.82 | 5×10⁻⁴ |
| DLPFC | 4916 | 1.80 | 1.77 | +0.61 | 5×10⁻⁴ |
| Hippocampus | 7799 | 2.09 | 2.07 | +0.69 | 5×10⁻⁴ |

\*Tech-adj currently includes coverage + CpG variance (LINE/mappability join incomplete for ~half of VMRs; matched analysis still matches on length + umap where available).

Mean VMR proportion of significant CpGs rose vs M0 (caudate 0.38→0.42; DLPFC 0.25→0.30; hippocampus 0.28→0.33), so absolute predictability coefficients are slightly attenuated but remain strongly positive in all regions.

Direction is consistent in all three regions; adjustment and matching preserve the effect.

Outputs: `{region}/vmr_meqtl_burden.tsv.gz`, `burden_model_results.tsv`, `matched_analysis_results.tsv`.

## Historical success criteria (not current evidence)

1. Burden increases with continuous predictability — **pass** (all regions)
2. Survives adjustment — **pass** (minimal + coverage/variance)
3. Survives matched analysis — **pass** (all regions)
4. Direction consistent in ≥2 regions — **pass** (3/3)
5. Supported by ≥1 external resource — **pass** (Jaffe/DLPFC + Schulz/hippocampus; see `../03_external_meqtl_validation/`)
