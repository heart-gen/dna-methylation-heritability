# Cell-type LINE/L1 sensitivity — results

**Decision:** `keep_main_figure_with_cell_composition_row`

LINE/L1–repressive enrichment of high-predictability VMRs remains after adjustment for cell-composition–correlated methylation properties in caudate and ≥1 other region.

## Claim summary

| Analysis | Regions dir. consistent (cellPC) | Regions p&lt;0.05 | Caudate pass | Claim pass |
|---|---:|---:|---|---|
| LINE_L1_enrichment_vs_predictability | 3 | 2 | True | True |
| H3K9me3_enrichment_vs_predictability | 3 | 3 | True | True |
| quiescent_chromatin_enrichment_vs_predictability | 3 | 3 | True | True |
| predictability_meqtl_burden_association | 3 | 3 | True | True |

## Key outputs

- `/projects/b1213/users/kynon/projects/dna-methylation-heritability/meqtl-validation/11_celltype_compartment_sensitivity/_m/enrichment_celltype_models.tsv`
- `/projects/b1213/users/kynon/projects/dna-methylation-heritability/meqtl-validation/11_celltype_compartment_sensitivity/_m/sample_level_celltype_checks.tsv`
- `/projects/b1213/users/kynon/projects/dna-methylation-heritability/meqtl-validation/11_celltype_compartment_sensitivity/_m/celltype_robustness_table.tsv`
- Phase 6 table updated: `/projects/b1213/users/kynon/projects/dna-methylation-heritability/meqtl-validation/07_repeat_mappability_sensitivity/_m/consolidated_robustness_table.tsv`

## Interpretation rules

- Genomic LINE/L1 overlap is unchanged by deconvolution; adjustment targets methylation variance explained by cell composition as a confounder of the predictability ranking.
- Significant sample-level Oligo ~ LINE-VMR methylation indicates bulk cellularity tracks repeat-rich methylation, but does not by itself refute sequence/compartment enrichment if cellPC-adjusted ORs remain &gt;1.
