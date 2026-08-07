# Cell-composition estimates

This module contains two deliberately separate cell-deconvolution tracks:

- `music-proportions-*.tsv`: existing RNA-seq MuSiC estimates. These are an
  orthogonal validation source and are not overwritten.
- `dnam-scmd-proportions-*.tsv`: WGBS estimates at the pinned scMD Lee and Tian
  850K single-cell DNAm signature CpGs (450K is a coordinate fallback).

## Analysis intent

We will test whether DNAm-reference cell fractions from each donor-region WGBS
sample agree with broad, independently estimated RNA cell fractions, supporting
their use as sensitivity-only covariates for CpG cis-meQTL mapping.

## Strategy memo

### 1. Biological Claim

DNAm-derived fractions estimate bulk-brain cellular mixture; they are not
direct histological cell counts.

### 2. Primary Hypothesis

DNAm and RNA total-neuronal fractions are positively correlated within each
region, with rho >= 0.30 and neuron-test FDR <= 0.05 required for downstream
integration.

### 3. Statistical Unit and Design

The unit is one adult donor-region sample. Caudate, DLPFC, and hippocampus are
processed separately. AA and EA samples are estimated; only the prespecified AA
subset enters the M6d meQTL sensitivity model.

### 4. Primary Comparison

DNAm total neurons (excitatory + inhibitory) are compared with RNA total
neurons (excitatory + inhibitory, or D1-SPN + D2-SPN + inhibitory in caudate).

### 5. Covariates and Exclusions

Markers require coverage >=5 in >=80% of eligible samples. Samples require
>=80% of retained marker coordinates. Ambiguous liftover, blacklist, common
C/T-SNP, and duplicated reference coordinates are excluded.

### 6. Expected Effect Direction

Broad neuronal and matched glial estimates should correlate positively across
modalities; absolute proportions need not agree.

### 7. Minimal Primary Analysis

Download and checksum the pinned scMD reference, select reference markers,
lift hg19 markers to hg38, extract HDF5-backed WGBS beta values, deconvolve,
validate against MuSiC, and write QC artifacts.

### 8. Sensitivity Analyses

Lee and Tian references and component algorithms are retained separately.
Technical correlations with marker coverage, pH, and PMI are reported.
For regions failing the primary 850K-signature gate, the coordinate-keyed scMD
WGBS references are run as a prespecified platform/reference sensitivity under
`_m/dnam_scmd_wgbs_reference/`; they cannot replace the primary result merely
because their concordance is larger.

### 9. Orthogonal Validation

Existing region-matched RNA MuSiC estimates provide the independent assay.

### 10. Main Figure-Worthy Result

`dnam_scmd_qc_figure.{pdf,png}` contains fraction distributions, neuronal
cross-modality concordance, a cell-class correlation heatmap, and marker
coverage diagnostics. It is intended as a supplementary QC figure.

### 11. Reviewer Objections and Responses

Reference mismatch is addressed by region-specific RNA validation; genetic and
mapping artifacts by marker exclusions; same-assay circularity by retaining
M3a as primary and using M6d only as sensitivity.

### 12. Implementation Notes

Run from `inputs/cell_proportions/_m`:

```bash
conda env create --prefix ./conda_env \
  -f ../../../conda_environments/cell_deconvolution.yml
sbatch ../_h/step_3_reference.sh
sbatch --dependency=afterok:<reference_job> ../_h/step_4_extract.sh
sbatch --dependency=afterok:<extract_job> ../_h/step_5_deconvolve.sh
sbatch --dependency=afterok:<deconvolution_job> ../_h/step_6_validate.sh
```

Equivalently, `bash ../_h/submit_dnam_scmd.sh` submits the dependency-linked
chain and writes `dnam_scmd_submission_manifest.tsv` in `_m`.

`step_3_reference.sh` installs the locked EnsDeconv, MIND, and scMD revisions
into `reference/R_libs`. Set `INSTALL_SCMD=0` only to exercise the explicitly
recorded reference-ensemble fallback.

The runner prefers `scMD::scMD()` when the exact upstream package is installed.
Otherwise it records and uses `scmd_reference_ensemble_v1`, combining RPC,
constrained projection, and simplex-constrained least squares across the same
Lee and Tian 850K signatures. No fallback is silent.

After a region passes the integration gate, submit
`meqtl-validation/01_cpg_meqtl_mapping/_h/step_m6d.sh`, followed by
`meqtl-validation/02_vmr_meqtl_burden/_h/step_m6d_burden.sh`. These jobs write
only to M6d sensitivity directories. They never replace locked M3a outputs.

The coordinate-keyed WGBS sensitivity is reproducible with:

```bash
ref=$(sbatch --parsable ../_h/step_wgbs_reference.sh)
run=$(sbatch --parsable --dependency=afterok:${ref%%;*} ../_h/step_wgbs_extract_deconvolve.sh)
sbatch --dependency=afterok:${run%%;*} ../_h/step_wgbs_validate.sh
```
