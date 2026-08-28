# DNA Methylation Heritability in the Human Brain

Analysis code for the manuscript:

> **Local SNP-explained methylation variation reveals genetically anchored and
> exposure-associated methylation architecture in the human brain**
>
> Alexis Bennett, Elisa Kain Johnson, Nia N. Terry, Jalil Hemphill,
> Kynon J.M. Benjamin†
>
> *bioRxiv* (2026). DOI: [10.64898/2026.06.05.730443](https://doi.org/10.64898/2026.06.05.730443).
>
> † Corresponding author: kynon.benjamin@northwestern.edu

---

## Supplementary Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20547606.svg)](https://doi.org/10.5281/zenodo.20547606)

Supplementary data generated in this study are available at
https://doi.org/10.5281/zenodo.20547606.

---

## Overview

This repository contains the full analysis pipeline used to partition
variably methylated regions (VMRs) by SNP-explained variation in 168 admixed
Black American adults across three brain regions: the caudate nucleus,
dorsolateral prefrontal cortex (DLPFC), and hippocampus. Analyses distinguish
genetically anchored methylation (concentrated in repressive chromatin and
repeat-rich regions) from exposure-associated methylation (enriched in
gene-proximal, accessible chromatin).

---

## Repository Structure

### Active revision (v2)

The AJHG revision is organized as numbered modules that run in dependency order.
See `AGENTS.md` for the governing rules and `MIGRATION_MANIFEST.tsv` for how each
legacy directory maps onto them.

| Directory | Description |
|---|---|
| `00_shared/` | Shared library: config, donor identity/alignment, chromosome ordering, run provenance |
| `01_vmr_catalog/` | Corrected VMR discovery and per-VMR methylation phenotypes |
| `02_local_genetic_variance/` | `h2_en_calibrated` — primary quantitative endpoint |
| `03_local_snp_prediction/` | Held-out local SNP prediction (secondary endpoint) |
| `04_repeat_repressive_architecture/` | Repeat-rich and repressive compartments (primary biology) |
| `05_cpg_meqtl_burden/` | CpG cis-meQTL burden gradient |
| `06_partitioned_heritability/` | S-LDSC partitioned heritability on the continuous local SNP contribution score |
| `07_transcription_splicing_coupling/` | Expression and splicing coupling |
| `08_region_donor_generalization/` | Cross-region and donor-group generalization |
| `09_schizophrenia_risk_application/` | Schizophrenia-risk application |
| `10_integrated_manuscript_outputs/` | Manuscript tables, figures, and number registry |
| `config/` | Shared configuration for the above |

Modules 01 and 02 have accepted runs. Modules 03–05 are implemented and
smoke-verified but have no accepted run. Modules 06–09 are scaffolded only.
Module 10 has Figures 1–2 implemented. Each module is gated on its upstream
module recording a passing acceptance gate.

### Legacy directories

Retained for old-versus-new comparison during the revision, and retired only once
`MIGRATION_MANIFEST.tsv` records a validated v2 replacement. **Results in these
trees are not valid for scientific use** — see `writing-notes/PIPELINE_AUDIT.md`,
in particular defect V1 (donor row misalignment invalidating every VMR set) and
E1 (`r_squared_cv` is an in-sample fit, not prediction accuracy).

| Directory | Description |
|---|---|
| `vmr-analysis/` | Legacy VMR identification → `01_vmr_catalog/` |
| `calibrated-simulation-analysis/` | Calibrated variance estimator → `02_local_genetic_variance/` |
| `local-snp-prediction/` | Legacy elastic-net prediction → `03_local_snp_prediction/` |
| `meqtl-validation/` | CpG meQTL mapping, burden, repeat/cell sensitivities, Phase 7 → modules 04–09 |
| `environmental-analysis/` | Environmental proxy associations (supplemental at most) |
| `simulation-analysis/` | Validation simulations and method comparisons |
| `sensitivity-analysis/` | Stacked/Venn/Sankey figures (withdrawn; depend on `r_squared_cv > 0.75`) |
| `qc_analysis/` | Quality control and replication cohort analysis |
| `inputs/` | Reference files and input data (not distributed; see Data Availability) |

---

## Data Availability

Raw genotype and DNA methylation data are available from dbGaP under
accession [phs000979.v3.p2](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000979.v3.p2).

Supplementary processed data (VMR calls, heritability estimates, and
summary statistics) are available on Zenodo:
https://doi.org/10.5281/zenodo.20547606.

---

## Software Requirements

| Tool | Use |
|---|---|
| R (≥4.4) | Data processing, VMR identification, visualization |
| Python (≥3.10) | Supporting scripts and data wrangling |
| [PLINK2](https://www.cog-genomics.org/plink/2.0/) | Genotype extraction and LD-based filtering |
| [GENBoostGPU](https://github.com/heart-gen/GENBoostGPU) | GPU-accelerated elastic-net SNP heritability estimation |
| [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) | GREML-based heritability comparison |
| Conda | Environment management (`epigenomics` env for R; `genomics` env for liftover) |

Pipeline steps are designed for SLURM-based HPC systems. Submission scripts
(`.sh`) are located in `_h/` subdirectories within each analysis module.

---

## Citation

If you use this code or data, please cite:

```
Bennett A, Johnson EK, Terry NN, Hemphill J, and Kynon JM Benjamin.
Local SNP-explained methylation variation reveals genetically anchored and
exposure-associated methylation architecture in the human brain.
bioRxiv (2026). DOI: 10.64898/2026.06.05.730443.
```
