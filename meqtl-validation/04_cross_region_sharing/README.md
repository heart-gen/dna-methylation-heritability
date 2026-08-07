# Phase 4a: Cross-region meQTL sharing

Test sharing and concordance of CpG-level cis-meQTL evidence and VMR meQTL support across caudate, DLPFC, and hippocampus (AA discovery).

## Inputs

- Phase 1 lead tables: `01_cpg_meqtl_mapping/{region}/_m/tensorqtl/qc/lead_snp_per_cpg.tsv.gz` (**M3a**)
- Phase 2 burden: `02_vmr_meqtl_burden/_m/{region}/vmr_meqtl_burden.tsv.gz` (**M3a**)
- Sample inclusion lists for shared-donor counts

## Analyses

1. VMR meQTL-support Jaccard / replication fractions (pairwise + all3)
2. Logistic: shared support ~ mean predictability
3. CpG lead-effect direction concordance and Pearson correlations (slope and z = slope/SE)
4. Discovery rates with √N normalization + equal-CpG-set contrasts
5. AA burden-gradient coefficients by region
6. Shared-donor counts (92 donors in all 3 regions)
7. **Experiment 3:** N-matched caudate lead-SNP retention downsample + shared-donor G×region

## Run

```bash
cd meqtl-validation/04_cross_region_sharing/_m
mkdir -p logs
sbatch ../_h/step_1.sh
# Experiment 3:
sbatch ../_h/step_3_downsample.sh   # N-matched lead-SNP retention (30 reps)
sbatch ../_h/step_4_gxregion.sh     # shared-donor G×region architecture screen
# after donor-group step_1/2 as well:
sbatch --dependency=afterok:<cross>:<down>:<gx>:<donor> ../_h/step_2.sh
```

Outputs for Experiment 3:
- `_m/caudate_downsample/` — design, replicate results, claim snapshot
- `_m/gxregion/` — pair results, claim snapshot
- Updates `pending_analyses.tsv` to `done` when complete

Method note: downsample uses **lead-SNP retention** under M3a residualization (same design as Phase 7 Tier A), not full cis remapping of each replicate.

## Key results (2026-08-01; M3a calls)

| Contrast | Shared CpGs | Both sig | Direction concordance | Pearson z (either sig) |
|---|---:|---:|---:|---:|
| Caudate–DLPFC | 69,030 | 25,345 | 0.912 | 0.894 |
| Caudate–Hip | 71,892 | 27,936 | 0.905 | 0.885 |
| DLPFC–Hip | 106,366 | 33,769 | 0.924 | 0.923 |

VMR support Jaccard is modest (~0.28–0.30 pairwise; 260 VMRs supported in all 3), but higher predictability predicts shared support (logistic OR ≈ 1.81–2.02, all pairwise p ≪ 0.05).

Prior M0-era outputs archived at `_m/M0_archive_20260730/`.  
Claim summary: `_m/phase4_claim_summary.tsv`.
