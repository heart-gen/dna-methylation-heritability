# Phase 3: External public brain meQTL validation

## Question

Are CpGs / VMRs with higher local genetic predictability more likely to have
independent published brain cis-meQTL evidence?

## Catalog

See `inputs/data_dictionary/_m/public_meqtl_resources.tsv` and
`inputs/data_dictionary/_h/README.md`.

| Resource | Role | Status |
|---|---|---|
| `jaffe_dlpfc_450k_meqtl` | Primary external (DLPFC array) | Harmonized + enrichment done |
| `schulz_hippocampus_array_meqtl` | Primary external (hippocampus array) | Harmonized + enrichment done |
| `brainseq_wgbs_meqtl` | **Not external** — cohort overlaps this discovery data (BrainSeq/LIBD) | Synapse deferred; even if obtained, do not count as independent Phase 3 validation |
| `brainseq_wgbs_meqtl_scz_subset` | Exploratory only (Nature SCZ-risk tables) | Harmonized; not genome-wide / not independent |

**Note:** BrainSeq full catalogs are same-/overlapping-cohort with this project's WGBS+genotype donors. Prefer Jaffe and Schulz for independent external support. See `_m/harmonized/brainseq_wgbs_meqtl.PENDING_SYNAPSE.txt`.


## Run

```bash
cd meqtl-validation/03_external_meqtl_validation/_m
mkdir -p logs
sbatch ../_h/step_1.sh          # init / download checklist
sbatch ../_h/step_harmonize.sh  # hg38 harmonize + VMR overlap
# After Phase 2 burden tables exist:
sbatch ../_h/step_2.sh          # enrichment models + summary
# Or locally:
#   conda activate /projects/p32505/opt/envs/genomics
#   python3 ../_h/04_run_and_summarize.py
```

Primary model: VMR any-external-support ~ continuous predictability (+ n CpGs / tech covariates).  
Preferred tissue pairings are `primary`; other region×resource tests are `secondary_cross_region`.

## Primary results (2026-07-30)

Adjusted-minimal binomial GLM (predictability + n_tested_cpgs):

| Resource | Region | Role | Support rate | Coef | p |
|---|---|---|---:|---:|---:|
| Jaffe 450K | DLPFC | primary | 0.24 | +0.15 | 3.5×10⁻⁶ |
| Schulz array | Hippocampus | primary | 0.06 | +0.30 | 2.4×10⁻¹⁷ |

Matched high vs low predictability (perm. p): Jaffe/DLPFC Δ=+0.11 (p≈0.004); Schulz/hippocampus Δ=+0.09 (p≈5×10⁻⁴).

**§7.5 criterion 5:** pass (≥1 external resource supports the gradient; both primary pairings pass).

Note: Jaffe × caudate cross-region secondary is negative — interpret tissue-matched primary tests for the claim; do not pool platforms.

Outputs:

- `{region}/external_support_model_{resource}.tsv`
- `{region}/external_matched_{resource}.tsv`
- `_m/external_support_models_all.tsv`
- `_m/external_support_primary_summary.tsv`
- `_m/phase3_criterion5_verdict.tsv`
