# Phase 6: Repeat / mappability technical validation

> **Repair-v2 order:** run `step_2_tech_joins.sh` after Phase 2 aggregation but
> before Phase 2 modeling. It runs in join-only mode, requires at least 50%
> reciprocal overlap for non-exact interval matches, records join provenance, and
> cannot update historical Phase 6 claim tables. Run `step_1.sh` after repaired
> Phase 2 models are complete.

Produce **one** consolidated robustness table (not scattered sensitivity figures).

## Question

Do LINE/L1, H3K9me3, and quiescent enrichments of high-predictability VMRs — and the
predictability→meQTL-burden association — survive mappability, segdup, SNP-proximity,
adjustment, and matching?

## Run

```bash
cd meqtl-validation/07_repeat_mappability_sensitivity/_m
mkdir -p logs
sbatch ../_h/step_0.sh   # annotation assets + VMR technical tables (if needed)
sbatch ../_h/step_2_tech_joins.sh   # join umap/LINE/SNP-prox/segdup onto Phase 2 burden
sbatch ../_h/step_1.sh   # consolidated robustness analyses
```

Primary script: `_h/03_run_robustness_analyses.py`. Tech joins: `_h/04_complete_tech_joins.py` → `_m/tech_join_completeness.tsv` (caudate 82%, hip 78%, DLPFC 48% EN∩burden ceiling).

## Sensitivities applied

| Filter | Rule |
|---|---|
| High mappability | Umap k24 mean ≥ 0.9 |
| Segdup excluded | `overlaps_segdup == 0` |
| SNP proximity excluded | no overlap with common SNP ±150 bp windows |
| Adjusted | length + num_snps (annotation) or coverage/variance (burden) |
| Matched | propensity-score high vs low predictability, 0.25-SD caliper, no replacement, paired randomization |

## Historical results (invalidated pending repair-v2 rerun; 2026-08-01 M3a burden arm)

Claim summary (`phase6_claim_summary.tsv`): direction-consistent in ≥2 regions.

| Analysis | Regions direction-consistent | Pass |
|---|---:|---|
| LINE/L1 ~ predictability | 2 / 3 | yes |
| H3K9me3 ~ predictability | 3 / 3 | yes |
| Quiescent ~ predictability | 3 / 3 | yes |
| Predictability → meQTL burden (M3a) | 3 / 3 | yes |

Prior outputs archived at `_m/M0_archive_20260730/`.

**Cell-composition sensitivity (done):** [`../11_celltype_compartment_sensitivity/`](../11_celltype_compartment_sensitivity/) — claim kept with cell-composition row; DLPFC LINE/L1 still fragile. Plan: [`CELLTYPE_LINE_L1_PLAN.md`](CELLTYPE_LINE_L1_PLAN.md).

Outputs:

- `_m/consolidated_robustness_table.tsv` — main table
- `_m/phase6_claim_summary.tsv`
- `_m/{region}/robustness_results.tsv`
- `_m/{region}/vmr_features_with_tech.tsv.gz`
