# Phase 5: Transcription and splicing integration

Reuse existing VMR–expression and VMR–PSI associations from `regulatory_context`. Do **not** rerun genome-wide transcriptome screens.

## Question

Are VMRs with CpG-level cis-meQTL support (and/or higher local genetic predictability) enriched for existing methylation–expression or methylation–PSI associations?

## Inputs

| Source | Path pattern |
|---|---|
| Phase 2 VMR meQTL burden | `02_vmr_meqtl_burden/_m/{region}/vmr_meqtl_burden.tsv.gz` |
| Expression (nearest gene, 250 kb) | `heritability/.../regulatory_context/_m/{region}/AA/expression/nearest_gene_window_250kb/architecture_model_input.tsv` |
| PSI (250 kb) | `.../AA/psi/window_250kb/architecture_model_input.tsv` |
| Expression ABC (secondary) | `.../AA/expression/abc/architecture_model_input.tsv` |

Join key: VMR genomic coordinates (`chrN_start_end` ↔ `N:start-end` via elastic-net summary `task_id`).

Outcome: `any_sig_fdr_05` from architecture tables.  
Exposure: any CpG with significant cis-meQTL (`n_cpgs_with_sig_meqtl > 0`); secondary predictors are continuous predictability and meQTL CpG proportion.

## Models

1. Fisher exact (one-sided greater): meQTL support × tx association  
2. Logistic (adjusted): `tx ~ meQTL_support + n_features_tested + vmr_length + min_distance + methylation_variance + num_snps`  
3. Logistic (adjusted): `tx ~ predictability + same covariates`  
4. Logistic (adjusted): `tx ~ meQTL_proportion + same covariates`  
5. Mann–Whitney: predictability in tx+ vs tx−  
6. Secondary: meQTL support × (expression ∩ PSI)

## Run

```bash
cd meqtl-validation/06_transcription_splicing_integration/_m
mkdir -p logs
sbatch ../_h/step_1.sh   # expression (array: 3 regions)
sbatch ../_h/step_2.sh   # PSI (+ optional ABC)
sbatch --dependency=afterok:<step1>:<step2> ../_h/step_3.sh
```

Scripts: `_h/01_meqtl_tx_enrichment.py`, `_h/02_summarize_phase5.py`.  
Env: `conda activate /projects/p32505/opt/envs/genomics`.

## Outputs

```text
_m/{region}/meqtl_x_{expression,psi,expression_abc}_enrichment.tsv
_m/{region}/vmr_meqtl_{modality}_joined.tsv.gz
_m/tx_enrichment_all.tsv
_m/tx_enrichment_primary.tsv
_m/tx_enrichment_both_modalities.tsv
_m/phase5_claim_summary.tsv
```

## Results (2026-08-01; M3a meQTL support)

| Region | Expression OR (p) | PSI OR (p) |
|---|---|---|
| Caudate | 1.88 (2.0e-8) | 4.80 (3.1e-25) |
| DLPFC | 12.78 (9.6e-20) | 1.89 (2.1e-3) |
| Hippocampus | 8.18 (6.6e-14) | 0.61 (0.89) |

Adjusted logistics for meQTL support agree in direction/significance for expression (all 3 regions) and PSI (caudate, DLPFC).

Prior M0-era outputs archived at `_m/M0_archive_20260730/`.

## Claim summary

- **Expression:** passes (3/3 regions significant positive).  
- **PSI:** passes (2/3 regions; hippocampus underpowered / not enriched).  
- Continuous predictability supports expression (DLPFC, hippocampus) and PSI (caudate, DLPFC).

Allowed manuscript framing: meQTL-supported VMRs show greater transcriptional and splicing coupling, strongest and most consistent for expression across regions and for PSI in caudate/DLPFC.
