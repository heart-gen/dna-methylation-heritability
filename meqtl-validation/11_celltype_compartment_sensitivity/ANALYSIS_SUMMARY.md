## Cell-type sensitivity for LINE/L1 and repressive-compartment enrichment

### Purpose

Test whether enrichment of high local-genetic-predictability VMRs in LINE/L1, H3K9me3, and quiescent chromatin—and the predictability→meQTL-burden association—survives adjustment for bulk cell-composition–correlated methylation properties in Black American (AA) caudate, DLPFC, and hippocampus. This is a Phase 6 sensitivity analysis supporting the repeat-rich / repressive-compartment architecture claim.

### Inputs

- AA meQTL inclusion covariates and CpG methylation beds: `meqtl-validation/01_cpg_meqtl_mapping/{region}/_m/prepared/`
- MuSiC RNA cell proportions and locked M5 `cellPC1–3`: `inputs/cell_proportions/_m/music-proportions-*.tsv`; `.../covariate_sensitivity/covariates_M5.txt`
- DNAm scMD proportions; caudate QC-gated M6d `dnamCellPC1–3`: `inputs/cell_proportions/_m/dnam-scmd-proportions-*.tsv`; `.../covariates_M6d.txt` (caudate only)
- Phase 6 VMR technical annotations: `07_repeat_mappability_sensitivity/_m/{region}/vmr_technical_annotations.tsv`
- Repeat / repressive annotation tables: `heritability/.../annotation/repeat_elements/_m/vmr_repeat_overlap_AA.tsv`, `.../repressive_chromatin/_m/vmr_repressive_overlap_AA.tsv`
- Phase 2 VMR meQTL burden: `02_vmr_meqtl_burden/_m/{region}/vmr_meqtl_burden.tsv.gz`

Sample overlap with AA meQTL donors: MuSiC 153/153 (caudate), 111/111 (DLPFC), 116/116 (hippocampus); DNAm scMD 151/153, 110/111, 116/116. DNAm integration QC gate passed for caudate only.

### Methods Text

We quantified, for each VMR with CpG-level methylation in the AA meQTL cohort, the mean methylation across CpGs within the VMR and the fraction of cross-sample methylation variance explained by bulk cell-composition covariates (`cellPC_r2` from OLS of VMR mean methylation on MuSiC `cellPC1–3`; caudate also `dnamCellPC_r2` from locked M6d DNAm composition PCs). Oligodendrocyte fraction (*Oligo*) and absolute Oligo–methylation correlation were recorded as secondary composition metrics.

Primary VMR-level models were binomial GLMs of binary genomic annotation (LINE/L1 overlap, H3K9me3 overlap, or quiescent chromatin) on continuous local genetic predictability, with covariates for VMR length and local SNP count (technical adjustment) and, in the primary cell-composition sensitivity, `cellPC_r2`. Parallel models adjusted for `oligo_r2`, absolute Oligo correlation, mean/variance of methylation, and—where available—`dnamCellPC_r2`. High-mappability subsets were evaluated as secondary filters. As a control outcome, we refit the predictability→meQTL-burden binomial model with the same cell-composition covariates. Predictors were z-scored within each model.

At the sample level, we tested associations between Oligo fraction or joint `cellPC1–3` and the mean methylation of LINE/L1 versus non-LINE VMRs (descriptive confounding check). Analyses were run separately by brain region. Claim success required direction-consistent cellPC-adjusted enrichment in ≥2 regions, cellPC-adjusted *P* < 0.05 in ≥2 regions, and a significant positive caudate cellPC-adjusted estimate.

Analyses were performed in Python using `statsmodels` GLMs/OLS and `scikit-learn` PCA (exploratory DNAm CLR-PCA in DLPFC/hippocampus only; primary DNAm inference restricted to caudate M6d). Deterministic seed `20260807` was used for exploratory DNAm PCA.

### Results Text

After joining annotation, technical, and cell-composition metrics, cellPC-adjusted models included 2,053 caudate, 520 DLPFC, and 1,689 hippocampal VMRs with complete data for LINE/L1 enrichment (530, 134, and 435 LINE/L1-overlapping VMRs, respectively).

**LINE/L1.** Continuous local genetic predictability remained positively associated with LINE/L1 overlap after `cellPC_r2` adjustment in caudate (OR = 1.55, *P* = 7.0 × 10⁻¹⁰) and hippocampus (OR = 1.20, *P* = 1.2 × 10⁻³), but not DLPFC (OR = 1.06, *P* = 0.54). Caudate DNAm `dnamCellPC_r2` adjustment agreed (OR = 1.48, *P* = 3.5 × 10⁻⁸). Relative to the unadjusted caudate OR (2.45), the cellPC-adjusted log-OR retained ~49% of the original signal (attenuation ratio 0.49). Under high-mappability restriction, caudate LINE/L1 enrichment remained significant (cellPC-adjusted OR = 1.60, *P* = 1.1 × 10⁻⁸), whereas DLPFC and hippocampus did not.

**H3K9me3 and quiescent chromatin.** CellPC-adjusted enrichment remained significant in all three regions for H3K9me3 (OR = 2.09 / 1.71 / 1.47; all *P* ≤ 1.9 × 10⁻⁴) and for quiescent chromatin (OR = 1.78 / 1.44 / 1.66; all *P* ≤ 2.3 × 10⁻⁴). Caudate DNAm adjustment was concordant for both marks.

**meQTL burden (control).** The predictability→meQTL-burden association remained strongly positive after `cellPC_r2` adjustment in all regions (binomial coefficients 3.51 / 2.20 / 2.34; all *P* ≈ 0).

**Sample-level confounding.** In caudate, Oligo fraction was associated with mean methylation of LINE/L1 VMRs (*β* = −0.057, *P* = 4.7 × 10⁻¹³, *R*² = 0.29) and non-LINE VMRs (*β* = −0.082, *P* = 2.1 × 10⁻¹⁷, *R*² = 0.38), indicating that bulk cellularity tracks methylation levels; the LINE-minus-nonLINE contrast was smaller but significant (*β* = 0.025, *P* = 9.6 × 10⁻⁵). Oligo associations with LINE-VMR mean methylation were not significant in DLPFC or hippocampus.

**Claim decision.** Prespecified criteria were met for LINE/L1 (pass in 2/3 regions), H3K9me3 (3/3), quiescent chromatin (3/3), and meQTL burden (3/3). The analysis decision was `keep_main_figure_with_cell_composition_row`: the LINE/L1–repressive enrichment claim is retained for the main manuscript with an explicit cell-composition sensitivity row, noting that DLPFC LINE/L1 remains fragile.

### Figure and Table Notes

- Potential main figure panel (Figure 5 robustness row): `meqtl-validation/07_repeat_mappability_sensitivity/_m/consolidated_robustness_table.tsv` columns `cellPC_adj_estimate`, `cellPC_adj_pvalue`, `dnamCellPC_adj_estimate`, `dnamCellPC_adj_pvalue`
  - Rationale: central technical defense of the repeat/repressive architecture claim
  - Key message: enrichment survives MuSiC cellPC_r2 adjustment in caudate (+ hippocampus for LINE/L1); DLPFC LINE/L1 does not
  - Required legend details: adjustment is for methylation variance explained by cell composition, not a change in genomic LINE/L1 overlap; DNAm PCs primary in caudate only

- Potential supplementary table: `meqtl-validation/11_celltype_compartment_sensitivity/_m/enrichment_celltype_models.tsv`
  - Rationale: full model grid (original, technical, cellPC, oligo, dnam, high-mappability)
  - Key columns: `region`, `analysis`, `model`, `or`, `pvalue`, `n`, `n_feature`

- Potential supplementary table: `.../_m/sample_level_celltype_checks.tsv`
  - Rationale: descriptive Oligo / cellPC confounding check
  - Key columns: `region`, `outcome`, `predictor`, `beta`, `pvalue`, `r2`

- Potential supplementary table: `.../_m/celltype_robustness_table.tsv` and `.../_m/celltype_claim_snapshot.tsv`
  - Rationale: manuscript-facing claim snapshot and attenuation ratios

No standalone figure PDF was generated by this workflow.

### Reproducibility Information

- Analysis directory: `meqtl-validation/11_celltype_compartment_sensitivity/`
- Primary scripts: `_h/01_audit_celltype_overlap.py`, `_h/02_build_vmr_cell_metrics.py`, `_h/03_run_celltype_sensitivity.py`, `_h/step_1.sh`
- Input files: listed under Inputs; sample overlap in `_m/celltype_sample_overlap.tsv`
- Output files: `_m/enrichment_celltype_models.tsv`, `_m/burden_celltype_models.tsv`, `_m/sample_level_celltype_checks.tsv`, `_m/celltype_robustness_table.tsv`, `_m/celltype_claim_summary.tsv`, `_m/celltype_claim_snapshot.tsv`, `_m/CELLTYPE_RESULTS.md`, `_m/{region}/vmr_mean_methylation.tsv.gz`, `_m/{region}/vmr_cell_metrics.tsv.gz`; Phase 6 table updated in place
- Log files inspected: `_m/logs/celltype_line.8779605.log`
- Execution command: `sbatch ../_h/step_1.sh` from `_m/` (SLURM job 8779605)
- Execution date: 2026-08-07 17:34:05–17:35:20 (from log)
- Git commit: `b1c99ed36` (short hash at summarization time; not recorded in the analysis log)
- Workflow manager: SLURM batch script (`step_1.sh`); no Snakemake/Nextflow
- Compute environment: Northwestern Quest HPC; partition `genomics`; account `b1042`; request 4 CPUs / 128G / 8 h
- Container or environment file: conda env `/projects/p32505/opt/envs/genomics` (activated in `step_1.sh`)
- R version: not used in this analysis
- Python version: 3.11.13 in `/projects/p32505/opt/envs/genomics` (queried at summarization; not printed in the SLURM log)
- Key package versions (same env query; not in SLURM log): `numpy` 2.2.6, `pandas` 2.3.3, `statsmodels` 0.14.5, `scikit-learn` 1.7.2
- Random seed: `20260807` (exploratory DNAm PCA in script `02_build_vmr_cell_metrics.py`)
- Missing reproducibility information: no per-job `pip freeze` in `celltype_line.8779605.log`; git commit not stamped in analysis outputs

### Limitations and Integration Notes

This analysis does not assign cell type of origin to individual VMRs and does not re-estimate elastic-net local genetic predictability after cell adjustment. Adjustment uses VMR-level methylation variance explained by deconvolution PCs as a proxy for composition-correlated methylation properties; genomic LINE/L1 overlap itself is unchanged by deconvolution. DNAm scMD PCs are primary only in caudate (integration QC gate); DLPFC/hippocampus DNAm CLR-PCA models are exploratory. DLPFC LINE/L1 enrichment was already fragile under high-mappability filters in Phase 6 and remains null after cellPC adjustment—manuscript text should treat caudate (+ hippocampus) as the supported LINE/L1 pattern and keep DLPFC as a caveat. Integrate with the Phase 6 consolidated robustness table (Figure 5) and avoid claiming that bulk cellularity is excluded as a contributor to methylation levels; the supported claim is that composition-correlated methylation properties do not abolish the predictability–compartment enrichment in the primary regions.
