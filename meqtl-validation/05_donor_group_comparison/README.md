# Phase 4b: Donor-group comparison

Primary discovery meQTL remains Black American (AA; locked **M3a**). White American (EA) stratified CpG meQTL uses:

- Genotypes: `inputs/genotypes/all_individuals/TOPMed_LIBD.*`
- CpG matrices: staff `vmr-analysis/all_individuals/{region}/_m/cpg/...`
- Covariates: M0 (`agedeath + sex + primarydx + snpPC1–5`); EA N too small for methPC lock
- Outputs: `01_cpg_meqtl_mapping/{region}/_m/{prepared,tensorqtl}/EA/`
- Burden: `02_vmr_meqtl_burden/_m/EA/{region}/`

## What this module does

1. Records EA meQTL readiness (`ea_meqtl_readiness.tsv`)
2. Tests AA vs EA local-predictability portability from elastic-net summaries
3. Lists cross-group concordant high-predictability VMRs + annotation / AA-meQTL enrichment
4. Compares AA vs EA burden-gradient coefficients
5. **Experiment 2 depth:** AA–EA CpG lead-effect concordance + MAF/cis-SNP-density matched discovery contrasts

Do **not** label AA–EA differences as ancestry-specific without formal interaction evidence and matched meQTL power.

## Run

```bash
# EA Phase 1 (does not overwrite AA primary):
cd meqtl-validation/01_cpg_meqtl_mapping/_m
sbatch --export=ALL,POPULATION=EA ../_h/step_1.sh
sbatch --dependency=afterok:<j1> --export=ALL,POPULATION=EA ../_h/step_2.sh
sbatch --dependency=afterok:<j1> --export=ALL,POPULATION=EA ../_h/step_3.sh
sbatch --dependency=afterok:<j2>:<j3> --export=ALL,POPULATION=EA ../_h/step_4.sh
sbatch --dependency=afterok:<j4> --export=ALL,POPULATION=EA ../_h/step_5.sh

# EA Phase 2:
cd meqtl-validation/02_vmr_meqtl_burden/_m
sbatch --export=ALL,POPULATION=EA ../_h/step_1.sh
sbatch --dependency=afterok:<agg> --export=ALL,POPULATION=EA ../_h/step_2.sh

# Refresh this module (portability + coefficients):
cd meqtl-validation/05_donor_group_comparison/_m
sbatch ../_h/step_1.sh
# Experiment 2 depth (concordance + MAF/LD matching):
sbatch ../_h/step_2.sh
```

EA-M3a caudate sensitivity (complete): `01_.../caudate/_m/tensorqtl/EA/M3a/` + `02_.../_m/EA_M3a/caudate/`.

## Key results (2026-08-01)

### Predictability portability

| Region | n shared | Pearson r | Spearman ρ | Jaccard high-pred (Q75) |
|---|---:|---:|---:|---:|
| Caudate | 11,466 | 0.87 | 0.51 | 0.57 |
| DLPFC | 9,593 | 0.42 | 0.12 | 0.22 |
| Hippocampus | 9,294 | 0.67 | 0.30 | 0.33 |

### EA stratified meQTL + burden

| Region | EA N | n FDR CpGs | Tech-adj burden coef | p |
|---|---:|---:|---:|---:|
| Caudate | 129 | 62,595 | 3.21 | ≪0.05 |
| DLPFC | 55 | 21,198 | 1.13 | ≪0.05 |
| Hippocampus | 60 | 29,372 | 1.29 | ≪0.05 |

Preferred interpretation: the aggregate predictability→meQTL-burden gradient is reproducible in both donor groups. Locus-level discovery counts differ with sample size (especially DLPFC/hippocampus EA).
