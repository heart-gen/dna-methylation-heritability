# LIBD / BrainSeq cis-eQTL mapping (AA)

Reusable genome-wide gene eQTL resource consumed by the Phase 7 Level 3 analysis.

This module is **not** schizophrenia-specific. SCZ risk-variant × gene tests live in
[`../08_schizophrenia_risk_application/`](../08_schizophrenia_risk_application/).

**Debug (low eGene count):** see [`EQTL_DEBUG_TODO.md`](EQTL_DEBUG_TODO.md).

## Locked recipe (caudate genes / standard)

| Item | Value |
|---|---|
| Cohort | `Age > 13`, `Race == AA`, `Dx ∈ {Control, SCZD}`, in `TOPMed_LIBD.AA` |
| Phenotype | edgeR `filterByExpr(~ Dx + Sex + Age)` → TMM `log2-CPM` |
| Covariates | `Sex + Dx + Age + snpPC1–3 + k` expression PCs (`k = num.sv`, capped at 50) |
| Cis window | ±500 kb (TSS) |
| MAF | 0.01 |
| Engine | TensorQTL cis permutation (`01_cpg_meqtl_mapping/_h/04_tensorqtl_map.py`) |
| Genotypes | `inputs/genotypes/TOPMed_LIBD.AA` |

## Run

Submit from `_m/`:

```bash
cd meqtl-validation/09_libd_eqtl_mapping/_m
mkdir -p logs
J1=$(sbatch --parsable --export=ALL,REGION=caudate ../_h/step_1.sh)   # prep + cov + BED
sbatch --dependency=afterok:${J1} --export=ALL,REGION=caudate ../_h/step_2.sh  # TensorQTL cis
```

## Outputs

```text
_m/<region>/genes/
├── prepared/          # expression BED, phenotypes, annotation
├── standard/          # covariates.txt, covariate_diagnostics.tsv
└── tensorqtl/         # *.cis_qtl.txt.gz
```

Caudate genes/standard (current lock): N=213, 33,206 genes, 20 expression PCs.
TensorQTL cis (2026-08-01): `libd_aa_caudate_genes_standard.cis_qtl.txt.gz`.

## Downstream consumers

- Phase 7 Level 3: `08_.../_h/17_libd_risk_variant_eqtl.py` (targeted index SNP × gene)
- Optional: transcription integration / locus panels

## Provenance

Initial caudate run was produced under `08_.../_m/level3/libd_eqtl/` and moved here
(`_m/REORGANIZED_FROM_08.tsv`). A compatibility symlink remains at
`08_.../_m/level3/libd_eqtl/caudate` → this module’s `_m/caudate`.

## Sensitivity: RPKM vs filterByExpr+CPM

Primary is under `genes/`. RPKM arm (mean RPKM > 0.2, `log2(RPKM+1)`, same covariates recipe)
writes to `genes_rpkm/` and does not overwrite primary:

```bash
cd meqtl-validation/09_libd_eqtl_mapping/_m
J1=$(sbatch --parsable --export=ALL,REGION=caudate ../_h/step_1_rpkm_sensitivity.sh)
J2=$(sbatch --parsable --dependency=afterok:${J1} --export=ALL,REGION=caudate ../_h/step_2_rpkm_sensitivity.sh)
sbatch --dependency=afterok:${J2} ../_h/step_3_compare_rpkm.sh
```

Comparison table: `_m/caudate/sensitivity_rpkm_vs_cpm.tsv`
