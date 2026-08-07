# Cell-type sensitivity for LINE/L1 / repressive-compartment claims

Implements [`../07_repeat_mappability_sensitivity/CELLTYPE_LINE_L1_PLAN.md`](../07_repeat_mappability_sensitivity/CELLTYPE_LINE_L1_PLAN.md).

## Question

Does enrichment of high local-genetic-predictability VMRs in LINE/L1, H3K9me3, and
quiescent chromatin survive adjustment for bulk cell-composition–correlated
methylation properties (MuSiC cellPCs; caudate DNAm scMD cellPCs)?

## Run

```bash
cd meqtl-validation/11_celltype_compartment_sensitivity/_m
mkdir -p logs
sbatch ../_h/step_1.sh
```

## Scripts

| Step | Script | Output |
|---|---|---|
| 1 | `01_audit_celltype_overlap.py` | `celltype_sample_overlap.tsv` |
| 2 | `02_build_vmr_cell_metrics.py` | `{region}/vmr_cell_metrics.tsv.gz`, meth matrices |
| 3 | `03_run_celltype_sensitivity.py` | enrichment + sample-level + claim snapshot; appends Phase 6 table |

## Primary tests

1. **VMR-level:** `annotation ~ predictability` with covariates for methylation variance
   explained by cellPCs (`cellPC_r2`) and Oligo correlation (`abs_oligo_corr`).
2. **Sample-level:** Oligo / cellPC association with mean methylation of LINE/L1 vs
   non-LINE VMRs.
3. **Burden sanity:** predictability → meQTL burden after adjusting for `cellPC_r2`.

DNAm scMD PCs are primary only for **caudate** (integration QC gate / M6d).
DLPFC/hippocampus use MuSiC; DNAm CLR-PCA is exploratory only (integration gate failed).

## Results (2026-08-07)

**Decision:** `keep_main_figure_with_cell_composition_row`

| Claim | cellPC-adj (Caud / DLPFC / Hip) | Pass |
|---|---|---|
| LINE/L1 | OR 1.55*** / 1.06 ns / 1.20** | Yes (2/3) |
| H3K9me3 | OR 2.09 / 1.71 / 1.47 (all sig) | Yes 3/3 |
| Quiescent | OR 1.78 / 1.44 / 1.66 (all sig) | Yes 3/3 |
| Burden | coef 3.51 / 2.20 / 2.34 (all sig) | Yes 3/3 |

See `_m/CELLTYPE_RESULTS.md` and `_m/celltype_claim_snapshot.tsv`.
