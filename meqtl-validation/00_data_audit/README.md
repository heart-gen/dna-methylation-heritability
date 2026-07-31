# Data Audit

## Analysis Summary

### Purpose

This primary Phase 0 audit asked whether the follow-up pilot has the required methylation, genotype, VMR, local genetic predictability, phenotype, covariate, cell-composition, expression, and PSI inputs needed before phenotype-wide implementation.

### Inputs

- Project root: `/projects/b1213/users/kynon/projects/dna-methylation-heritability`
- Staff project root: `/projects/b1213/users/alexis/projects/dna-methylation-heritability`
- Input configuration: `meqtl-validation/00_data_audit/_h/analysis_inputs.json`
- Phenotype and covariate table: `sample_summary/_m/phenotype_data.tsv`
- VMR BED files: `vmr-analysis/{caudate,dlpfc,hippocampus}/_m/vmr.bed`
- Local genetic predictability summaries: `heritability/elastic_net_model/all_individuals/{region}/_m/{region}_summary_elastic-net_AA.tsv`
- Cell-composition estimates: `inputs/cell_proportions/_m/music-proportions-*.tsv`
- Genotype symlinks: `inputs/genotypes/TOPMed_LIBD.AA.{pgen,psam,pvar,eigenvec,eigenval}`
- Molecular-link coverage files: `heritability/elastic_net_model/all_individuals/tissue_comparison/regulatory_context/_m/**/sample_feature_coverage.tsv`
- Staff CpG methylation matrices: `/projects/b1213/users/alexis/projects/dna-methylation-heritability/vmr-analysis/{region}/_m/cpg/chr_*/cpg_meth.phen`

### Methods Text

We performed a deterministic file and sample audit before phenotype-wide modeling. The audit enumerated required input files in the project and staff repositories, counted predefined VMRs and local genetic predictability summaries, summarized phenotype sample counts by region, donor group, and diagnosis, quantified missingness for primary and sensitivity covariates, recorded repeated-donor regional coverage, checked genotype symlink availability, inventoried CpG methylation matrices without scanning full wide matrices, and summarized existing VMR-expression and VMR-PSI link coverage.

No phenotype association models, meQTL models, or multiple-testing procedures were run in this audit.

### Results Text

The AA discovery phenotype table contained the expected region-specific sample counts: caudate 153 samples (88 control, 65 schizophrenia), DLPFC 111 samples (68 control, 43 schizophrenia), and hippocampus 116 samples (68 control, 48 schizophrenia). The expanded EA cohort was also represented: caudate 129 samples, DLPFC 55 samples, and hippocampus 60 samples.

The primary covariate set for schizophrenia and aging screening (`agedeath`, `sex`, `primarydx`, and `snpPC1`-`snpPC5`) was complete in all audited AA and EA region groups. Across all donor groups in `sample_summary/_m/phenotype_data.tsv`, 137 donors had all three regions represented, 43 had two regions, and 127 had one region.

Predefined VMR files were present in both repositories, with 12,001 caudate VMRs, 10,372 DLPFC VMRs, and 10,216 hippocampus VMRs. AA local genetic predictability summaries were present for 11,467 caudate VMRs, 9,594 DLPFC VMRs, and 9,295 hippocampus VMRs.

The local repo did not contain expanded per-chromosome CpG methylation matrices under `vmr-analysis/*/_m/cpg/chr_*`, but Alexis's repo did. In Alexis's repo, all three regions had 24 chromosome directories with both `cpg_meth.phen` and `res_cpg_meth.phen` files. Example chromosome 1 non-residualized matrices contained 1,994,366 columns for caudate, 1,810,615 columns for DLPFC, and 1,829,023 columns for hippocampus, including `FID` and `IID`.

AA genotype symlinks were available in the project repo. `TOPMed_LIBD.AA.psam` contained 525 samples and three columns; `TOPMed_LIBD.AA.pgen` and `TOPMed_LIBD.AA.pvar` were present via symlink to shared resources.

Existing VMR-expression and VMR-PSI link coverage files were available for all three regions. In the AA-only context, linked nearest-gene expression features numbered 4,050 in caudate, 2,915 in DLPFC, and 3,067 in hippocampus; linked PSI features numbered 426,144 in caudate, 350,632 in DLPFC, and 373,097 in hippocampus.

### Figure and Table Notes

- Potential supplementary table: `meqtl-validation/00_data_audit/_m/phenotype_counts.tsv`
  - Rationale: documents exact region, diagnosis, and donor-group sample sizes before modeling.

- Potential supplementary table: `meqtl-validation/00_data_audit/_m/primary_model_complete_cases.tsv`
  - Rationale: documents that the primary covariate set is complete for all audited region and donor-group strata.

- Potential supplementary table: `meqtl-validation/00_data_audit/_m/cpg_methylation_inventory.tsv`
  - Rationale: identifies that staff-repo CpG methylation matrices are needed for CpG-level meQTL mapping.

- Potential supplementary table: `meqtl-validation/00_data_audit/_m/molecular_link_coverage.tsv`
  - Rationale: records available expression and PSI link sets for later phenotype-specific convergence analyses.

No figure is recommended from the audit itself.

### Reproducibility Information

- Analysis directory: `meqtl-validation/data_audit`
- Primary script: `meqtl-validation/00_data_audit/_h/00_data_audit.py`
- Configuration files: `meqtl-validation/00_data_audit/_h/analysis_inputs.json`, `meqtl-validation/00_data_audit/_h/model_formulas.json`, `meqtl-validation/00_data_audit/_h/multiple_testing_families.tsv`
- Output files: `meqtl-validation/00_data_audit/_m/*.tsv`
- Log files inspected: none; this audit writes `meqtl-validation/00_data_audit/_m/audit_reproducibility.tsv`
- Execution command: `python3 meqtl-validation/00_data_audit/_h/00_data_audit.py`
- Execution date: recorded in `meqtl-validation/00_data_audit/_m/audit_reproducibility.tsv`
- Git commit: not recorded because `git status` failed under the sandbox when Git LFS tried to write inside `.git/lfs/tmp`
- Workflow manager: none
- Compute environment: local shell on the shared project filesystem
- Container or environment file: not used for this audit
- Python version: recorded in `meqtl-validation/00_data_audit/_m/audit_reproducibility.tsv`
- Random seed: not used; deterministic file audit
- Missing reproducibility information: no Git commit hash, no formal workflow execution log

### Limitations and Integration Notes

This audit confirms file availability and high-level sample readiness but does not validate genome-build consistency, CpG coordinate sorting, genotype-methylation ID concordance, or meQTL model calibration. The next implementation step should be a Phase 1 CpG-level meQTL preflight that joins staff-repo CpG methylation matrices to AA genotype sample IDs and predefined VMR intervals before launching region-specific cis-meQTL jobs.

The audit supports proceeding to Phase 1 CpG-level meQTL validation. It does not yet support choosing schizophrenia, aging, strengthened architecture, or stopping.
