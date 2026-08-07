# Phase 2: VMR-level aggregation of CpG meQTL evidence

## Question

Does CpG meQTL burden increase with continuous VMR local genetic predictability?

## Run

```bash
cd meqtl-validation/02_vmr_meqtl_burden/_m
mkdir -p logs
sbatch ../_h/step_1.sh   # aggregate (qval≤0.05 significance; joins coverage + technical annotations)
sbatch ../_h/step_2.sh   # burden models + matching
```

Primary outcome: proportion of CpGs with `qval ≤ 0.05`.  
Primary predictor: continuous VMR local genetic predictability (`h2_unscaled`).

**Phase 1 lock:** primary calls are **M3a** (`../01_cpg_meqtl_mapping/PHASE1_LOCK_DECISION.md`).  
Prior M0 Phase 2 outputs archived at `_m/M0_archive_20260730/`.

## Primary results (2026-08-01; M3a calls)

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

## Success criteria

1. Burden increases with continuous predictability — **pass** (all regions)
2. Survives adjustment — **pass** (minimal + coverage/variance)
3. Survives matched analysis — **pass** (all regions)
4. Direction consistent in ≥2 regions — **pass** (3/3)
5. Supported by ≥1 external resource — **pass** (Jaffe/DLPFC + Schulz/hippocampus; see `../03_external_meqtl_validation/`)
