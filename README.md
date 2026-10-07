# DNA Methylation Heritability in the Human Brain

Analysis code for the manuscript. The citation below is the v1 preprint; the
revision reframes it (see Overview) and its title will change:

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

This repository contains the analysis pipeline for a study of how local genetic
variation organizes interindividual DNA methylation variability across the
caudate nucleus, dorsolateral prefrontal cortex (DLPFC) and hippocampus, using
whole-genome bisulfite sequencing in admixed Black American adults. Local
genetic control is treated as a continuous, within-region relative rank of local
SNP contribution to variably methylated region (VMR) methylation; it is related
to repeat-rich and repressive chromatin, CpG meQTL burden, transcription and
splicing coupling, and trait-associated GWAS loci. Low local genetic control is
not read as evidence of environmental determination.

---

## Repository Structure

The analysis is organized as numbered modules that run in dependency order
(AGENTS.md 6). Each module records its accepted runs under **Accepted runs** in
its own README; that table, not this one, is the record of what may be cited.

| Directory | Description |
|---|---|
| `00_shared/` | Shared library: config, donor identity/alignment, chromosome ordering, run provenance, acceptance gates |
| `01_vmr_catalog/` | Corrected VMR discovery and per-VMR methylation phenotypes |
| `01b_estimation_cells/` | Donor-group and donor-count estimation cells on a pooled-discovery catalog |
| `02_local_genetic_variance/` | Relative local SNP contribution score (`local_snp_contribution_score`) -- the primary endpoint; absolute locus PVE is retired |
| `02b_greml_simulation_benchmark/` | GCTA-GREML recovery of absolute local h2 on simulated phenotypes; reads nothing from 02 |
| `02c_cis_greml_sensitivity/` | Conventional cis-GREML on the observed phenotypes against the 02 score: existence and ordering only, no per-VMR h2 |
| `03_local_snp_prediction/` | End-to-end out-of-fold local SNP prediction (secondary endpoint) |
| `04_repeat_repressive_architecture/` | Repeat-rich and repressive compartments (primary biology) |
| `05_cpg_meqtl_burden/` | CpG cis-meQTL burden gradient |
| `06_partitioned_heritability/` | S-LDSC on the continuous score, conditional on VMR membership |
| `07_transcription_splicing_coupling/` | Expression and splicing coupling |
| `08_region_donor_generalization/` | Tiered cross-region and donor-group generalization |
| `09_schizophrenia_risk_application/` | Schizophrenia-risk application and the GWAS negative-control collection |
| `09b_aging_application/` | Age-associated methylation differences along the axis |
| `10_environmental_exploratory/` | Exploratory exposure associations (supplement only) |
| `11_integrated_manuscript_outputs/` | Figures, tables, number registry, Methods/Results summaries |
| `config/` | Shared, PI-locked configuration |
| `inputs/` | Reference files and input data (not distributed; see Data Availability) |
| `supplementary_data/` | The deposition list: every Supplementary Data item and the accepted run that produced it |

### Retired v1 trees

The v1 analysis directories (`vmr-analysis/`, `calibrated-simulation-analysis/`,
`local-snp-prediction/`, `meqtl-validation/`, `environmental-analysis/`,
`simulation-analysis/`, `sensitivity-analysis/`, `qc_analysis/`,
`sample_summary/`, and the untracked `simulation-analysis.bak/`) were removed
from the working tree on 2026-10-06, after every row of `MIGRATION_MANIFEST.tsv`
was closed: a validated v2 replacement, a withdrawal, or a recorded PI decision
not to migrate. Their tracked files are recoverable from the annotated tag
**`v1-legacy-final`** (`git show v1-legacy-final:<path>`). The whole trees,
untracked outputs included, were moved to
`/projects/b1213/users/kynon/archive/dna-methylation-heritability-v1-20261003/`;
`legacy_v1_archive_inventory.tsv` summarizes them by subtree, and the full
per-file inventory (path, bytes, SHA-256, tracked flag) is in the archive's
`_inventory/`. **Results in those trees are
not valid for scientific use** -- see `writing-notes/PIPELINE_AUDIT.md`, in
particular defects V1 (donor row misalignment invalidating every VMR set) and E1
(`r_squared_cv` is an in-sample fit, not prediction accuracy).

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
