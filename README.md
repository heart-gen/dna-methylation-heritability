# DNA Methylation Heritability in the Human Brain

Analysis code for the manuscript:

> **Local SNP-explained methylation variation reveals genetically anchored and
> exposure-associated methylation architecture in the human brain**
>
> Alexis Bennett, Elisa Kain Johnson, Nia N. Terry, Jalil Hemphill,
> Kynon J.M. Benjamin†
>
> *bioRxiv* (2026). DOI: *pending screening*
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

| Directory | Description |
|---|---|
| `heritability/` | Elastic-net SNP heritability pipeline (VMR identification, GENBoostGPU estimation) |
| `vmr-analysis/` | Variably methylated region characterization and annotation |
| `environmental-analysis/` | Environmental proxy association analyses |
| `simulation-analysis/` | Validation simulations and method comparisons |
| `sensitivity-analysis/` | Sensitivity analyses (all individuals and Black American subgroup) |
| `covar-analysis/` | Covariate selection and testing |
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
| R (≥4.1) | Data processing, VMR identification, visualization |
| Python (≥3.9) | Supporting scripts and data wrangling |
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
Bennett A, Johnson EK, Terry NN, Hemphill J, Benjamin KJM.
Local SNP-explained methylation variation reveals genetically anchored and
exposure-associated methylation architecture in the human brain.
bioRxiv (2026). DOI: pending
```
