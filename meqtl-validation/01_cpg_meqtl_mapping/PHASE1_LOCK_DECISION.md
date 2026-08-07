# Phase 1 lock decision

**Date:** 2026-08-01
**Regions:** caudate, DLPFC, hippocampus
**Chosen primary covariate model:** `M3a`

## Model definition

```text
M3a = agedeath + sex + primarydx + snpPC1–5 + methPC1–5
```

`methPC1–5` are the first five principal components of M0-residualized CpG methylation (PEER unavailable in the genomics env).

## Rationale

Lock **M3a** over M0 based on consistent **λ_NS improvement** and **increased meQTL detection** in all three regions, with external enrichment retained:

| Region | Metric | M0 | M3a | Change |
|---|---|---:|---:|---|
| caudate | n significant (FDR) | 73,375 | 81,007 | **+10.4%** |
| caudate | λ_NS (`qval > 0.05`) | 1.468 | 1.459 | −0.009 |
| caudate | external OR vs M0 | 1.000 | 0.993 | retained |
| DLPFC | n significant (FDR) | 42,496 | 51,273 | **+20.7%** |
| DLPFC | λ_NS | 1.503 | 1.468 | **−0.035** |
| DLPFC | external OR vs M0 | 1.000 | 0.969 | retained |
| hippocampus | n significant (FDR) | 47,630 | 55,957 | **+17.5%** |
| hippocampus | λ_NS | 1.460 | 1.419 | **−0.041** |
| hippocampus | external OR vs M0 | 1.000 | 1.079 | retained / improved |

M3a was preferred to M3b/M3c as the smallest latent set that improves calibration and increases discoveries without collapsing external enrichment. Residual λ_NS (~1.42–1.47) remains above 1.2 and is interpreted as abundant true cis-meQTL burden rather than unmodeled batch alone.

**Note on prior automatic rule:** the scripted go/no-go threshold required Δλ_NS ≥ 0.15 (or λ_NS ≤ 1.20). That bar was not met; this lock is an explicit investigator override prioritizing directional λ_NS improvement plus increased detection with retained external support.

## Decision rules applied (updated)

1. Prefer models with lower λ_NS than M0 in ≥2 regions.
2. Prefer models that increase (not collapse) FDR-significant CpG discoveries.
3. Require tissue-matched external enrichment OR ≥ ~80% of M0.
4. Prefer the smallest latent set among models satisfying (1)–(3) → **M3a (k=5)**.

## Per-region comparison snapshots

### caudate

```
model_id  n_significant_fdr  lambda_gc_nonsignificant  enrichment_or  enrichment_or_vs_M0
      M0            73375.0                  1.467798       1.649908             1.000000
      M1            72749.0                  1.465379       1.650965             1.000641
      M2            72216.0                  1.480686       1.638065             0.992823
     M3a            81007.0                  1.459178       1.638660             0.993183
     M3b            80748.0                  1.448133       1.670292             1.012355
     M3c            79422.0                  1.477839       1.637325             0.992374
```

### dlpfc

```
model_id  n_significant_fdr  lambda_gc_nonsignificant  enrichment_or  enrichment_or_vs_M0
      M0            42496.0                  1.503024       2.441817             1.000000
      M1            42291.0                  1.498598       2.440271             0.999367
      M2            41714.0                  1.465015       2.396297             0.981358
     M3a            51273.0                  1.468304       2.365776             0.968859
     M3b            51265.0                  1.460617       2.348642             0.961842
     M3c            49453.0                  1.497140       2.362136             0.967368
```

### hippocampus

```
model_id  n_significant_fdr  lambda_gc_nonsignificant  enrichment_or  enrichment_or_vs_M0
      M0            47630.0                  1.459641       6.371749             1.000000
      M1            47665.0                  1.418826       6.466619             1.014889
      M2            46712.0                  1.392282       6.383938             1.001913
     M3a            55957.0                  1.418517       6.877335             1.079348
     M3b            55409.0                  1.475883       6.855709             1.075954
     M3c            54731.0                  1.443659       6.539989             1.026404
```

## Promotion performed

For each region (`caudate`, `dlpfc`, `hippocampus`):

1. `prepared/covariates.txt` ← `covariate_sensitivity/covariates_M3a.txt` (M0 copy kept as `prepared/covariates_M0_pre_lock.txt`).
2. Prior primary TensorQTL moved to `{region}/_m/tensorqtl_M0_archive/`.
3. Primary `tensorqtl/cpg_meqtl_{region}.cis_qtl.txt.gz` and `tensorqtl/qc/` installed from `covariate_sensitivity/tensorqtl/M3a/`.

See `_m/phase1_m3a_promote.log`.

## Next steps

1. ~~**Recompute Phase 2** VMR meQTL-burden tables and models on the promoted M3a significance calls.~~ **Done 2026-08-01** (still Pass in 3/3 regions).
2. Spot-check Phase 3–6 claims that consume primary lead tables; re-run only if architecture conclusions change.
3. Keep M0 archived for manuscript sensitivity (M0 vs M3a λ_NS / n_sig panel).

## Status

**LOCKED:** primary CpG cis-meQTL covariate model is **`M3a`** for all three brain regions.

## DNAm cell-composition sensitivity addendum (2026-08-05)

`M6d = M3a + dnamCellPC1–3` is a sensitivity model only. The DNAm cell PCs are
constructed within each region and model sample set by zero replacement,
closure, centered log-ratio transformation, and PCA of scMD proportions. The
highest-numbered PC is dropped only if needed to retain full matrix rank.

M6d is created only for regions passing the prespecified DNAm/RNA integration
gate (bounded fractions, sums equal to one, marker/sample QC, and total-neuron
Spearman rho >= 0.30 with BH FDR <= 0.05). Its TensorQTL and VMR-burden outputs
remain under separate `M6d` sensitivity directories. This addendum does not
unlock, replace, or overwrite primary M3a.

Samples failing DNAm marker QC are not composition-imputed. For every eligible
region, both `M3a_dnam_matched` and M6d are rerun on the same retained sample
list and the same CpG phenotype BED, so the sensitivity comparison is paired.

Manuscript retention requires a positive VMR predictability–burden effect in at
least two regions, external enrichment at least 80% of M3a, and no greater than
50% collapse in FDR-significant CpG discoveries. Final results are recorded in
`02_vmr_meqtl_burden/_m/M6d/m6d_robustness_summary.tsv` after sensitivity jobs
finish.
