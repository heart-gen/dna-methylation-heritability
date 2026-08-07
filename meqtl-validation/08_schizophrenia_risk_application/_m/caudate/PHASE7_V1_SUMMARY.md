# Phase 7 v1 summary — caudate AA (M3a)

**Date:** 2026-08-01  
**Jobs:** `8457333` (define+link), `8458058` (targeted meQTL + architecture + tx; sample-order fix)

## Results

| Metric | Value |
|---|---:|
| Index SNPs in genotype panel | 489 / 624 |
| Index SNPs tested (MAF≥0.05) | 377 |
| Risk-variant–CpG tests | 40,001 |
| FDR<0.05 pairs | 361 |
| Significant CpGs / VMRs / loci | 361 / 38 / 31 |
| Primary VMR links (Analysis 2) | 865 VMRs / 310 loci |
| Architecture tech-adj coef (predictability) | 0.526 (p=1.6×10⁻⁴) |
| Matched high−low Δ(any_sig) | 0.096 (perm p=0.003) |
| SCZ-meQTL VMRs with expression link | 4 (Fisher p=0.022) |
| SCZ-meQTL VMRs with PSI link | 4 (Fisher p=0.001) |

## Implementation note

First `step_2` run (`8457334`) returned 0 FDR hits due to a `PgenReader(select_samples=...)` sample-order bug (dosages stayed in genotype order while phenotypes used BED order). Fixed in `03_test_risk_variant_cpg_meqtl.py` before `8458058`.

## Deferred

- Analysis 5: regional specificity / caudate downsampling
- Analysis 7: diagnosis association at prioritized loci
- Formal colocalization
- LD-proxy R² expansion beyond index-SNP identity
