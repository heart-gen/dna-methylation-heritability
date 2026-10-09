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

The analysis is organized as numbered modules that run in dependency order. Each module records its accepted runs under **Accepted runs** in
its own README; that table, not this one, is the record of what may be cited.

The numbers give a run order, not a chain. The graph shows what each module
reads:

```mermaid
flowchart TD
    m01["01 VMR catalog"]
    m01b["01b estimation cells"]
    m02["02 local SNP contribution score"]
    m02b["02b GREML simulation benchmark"]
    m02c["02c cis-GREML sensitivity"]
    subgraph score["Read the 02 score"]
        m03["03 local SNP prediction"]
        m04["04 repeat and repressive architecture"]
        m05["05 CpG meQTL burden"]
        m06["06 partitioned heritability"]
        m07["07 transcription and splicing coupling"]
        m08["08 region and donor generalization"]
        m09["09 schizophrenia GWAS loci"]
        m09b["09b aging"]
        m10["10 environmental, exploratory"]
    end
    m11["11 manuscript figures and tables"]

    m01 --> m01b --> m02
    m01 --> m02
    m01 --> m02b
    m02 --> m02c
    m02 --> score
    m03 --> m04
    m04 --> m08
    m04 --> m09b
    m04 --> m10
    m05 --> m07
    m07 --> m08
    m07 --> m09b
    m06 --> m09
    m08 --> m09
    score --> m11
    m02b --> m11
    m02c --> m11

    classDef aside stroke-dasharray: 5 5
    class m02b,m02c,m10 aside
```

Each arrow is an upstream run ID recorded in an accepted run's `manifest.tsv`.
Two kinds of arrow are left out to keep the graph readable:

- 01 into the modules that also read 02. Every module except 06 reads the 01
  catalog directly.
- Arrows implied by a longer path: 03 also reads 01b; 08 also reads 01b, 03 and
  05; 09 also reads 04, 05 and 07.

Dashed boxes are the benchmark, sensitivity and exploratory modules. Update the
graph when a module gains or drops an upstream.

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
| `09_schizophrenia_gwas_loci/` | Schizophrenia GWAS loci on the local-genetic-control axis, and the GWAS negative-control collection |
| `09b_aging_application/` | Age-associated methylation differences along the axis |
| `10_environmental_exploratory/` | Exploratory exposure associations (supplement only) |
| `11_integrated_manuscript_outputs/` | Figures, tables, number registry, Methods/Results summaries |
| `config/` | Shared, PI-locked configuration |
| `inputs/` | Reference files and input data (not distributed; see Data Availability) |
| `supplementary_data/` | The deposition list: every Supplementary Data item and the accepted run that produced it |

### Retired v1 trees

The v1 analysis directories (`vmr-analysis/`, `calibrated-simulation-analysis/`,
`local-snp-prediction/`, `meqtl-validation/`, `environmental-analysis/` and
others) were retired on 2026-10-06 after each had a validated v2 replacement, a
withdrawal, or a recorded decision not to migrate. Their tracked files are
recoverable from the annotated tag **`v1-legacy-final`**
(`git show v1-legacy-final:<path>`). **Results in those trees are not valid for
scientific use**: the v1 VMR sets carry a donor-row misalignment, and the v1
`r_squared_cv` is an in-sample fit, not prediction accuracy.

---

## Data Availability

Raw genotype and DNA methylation data are available from dbGaP under
accession [phs000979.v3.p2](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000979.v3.p2).
No individual-level methylation, genotype or covariate data are distributed in
this repository or on Zenodo.

Results are distributed in three tiers, recorded file by file in
`supplementary_data/release_manifest.tsv` and built by
`supplementary_data/_h/build_release.py`:

| Tier | What | Where |
|---|---|---|
| git | Every file the manuscript is written from, up to 15 MB: figure source data, Supplementary Data tables, decision and gate tables, run manifests, the accepted figure run | this repository |
| Git LFS | The same class of file above 15 MB | this repository (Git LFS) |
| Zenodo | Reproducibility extras no manuscript text uses (per-run provenance, uncited per-VMR tables, annotations, LD scores), plus a mirror of the two tiers above; one zip per module | [10.5281/zenodo.20547606](https://doi.org/10.5281/zenodo.20547606) |

Bulk intermediates (per-task shards, full nominal meQTL output, coloc region
shards, copies of public GWAS) are regenerable from each module's `_h/` and are
not deposited.

---

## Software Requirements

| Tool | Use |
|---|---|
| R (≥4.4) | Data processing, VMR identification, visualization |
| Python (≥3.10) | Supporting scripts and data wrangling |
| [PLINK2](https://www.cog-genomics.org/plink/2.0/) | Genotype extraction and LD-based filtering |
| [glmnet](https://glmnet.stanford.edu/) | Elastic-net local SNP models (Modules 02, 03) |
| [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) | GREML benchmark and cis-GREML sensitivity (Modules 02b, 02c) |
| [TensorQTL](https://github.com/broadinstitute/tensorqtl) | CpG cis-meQTL mapping (Module 05) |
| [LDSC](https://github.com/bulik/ldsc) | Stratified LD score regression (Module 06) |
| [coloc](https://chr1swallace.github.io/coloc/) / [susieR](https://stephenslab.github.io/susieR/) | Colocalization (Module 09) |
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
