# DNAm cell-deconvolution QC decision

**Date:** 2026-08-05  
**Primary signatures:** scMD Lee and Tian 850K, package commit
`9b4f52f4bc22ab1e39266f61a80233597e6b52c1`  
**Genome harmonization:** EPIC hg19 annotation 0.6.0, one-to-one liftOver to hg38  
**Primary model status:** M3a remains locked and unchanged

## Executed workflow

Raw methylated and total counts were extracted at eligible signature CpGs from
the three HDF5-backed BSseq objects. Marker coordinates required coverage >=5
in >=80% of samples; samples required >=80% retained-marker coverage. Isolated
missing beta values were median-imputed by marker within region. Two caudate,
one DLPFC, and zero hippocampal samples failed and were excluded without
sample-level composition imputation.

The active legacy epigenomics R prefix cannot compile upstream scMD dependencies
because its Makeconf lacks `SHLIB_LIBADD`. This failure is written to
`_m/reference/scmd_install_status.tsv`. The completed estimates therefore use
the explicitly named `scmd_reference_ensemble_v1`: RPC, constrained projection,
and simplex-QP components across the same pinned Lee and Tian 850K signatures.
Component and reference-specific estimates are retained. No result is labeled
as an upstream `scMD::scMD()` run.

## Prespecified integration gate

| Region | Estimated / expected | Total-neuron Spearman rho | Neuron FDR | Pass |
|---|---:|---:|---:|:---:|
| Caudate | 280 / 282 | 0.720 | 2.05e-45 | Yes |
| DLPFC | 165 / 166 | -0.035 | 0.659 | No |
| Hippocampus | 176 / 176 | -0.103 | 0.265 | No |

All estimated fractions were bounded and summed to one. Lee- and Tian-specific
total-neuron estimates gave the same regional conclusions.

## Coordinate-keyed WGBS reference sensitivity

Because DLPFC and hippocampus failed with the primary 850K signatures, scMD's
coordinate-keyed Lee/Tian WGBS references were tested as a platform/reference
sensitivity using identical coverage, exclusion, imputation, and RNA-validation
rules. This was scientifically prespecified and was not selected based on its
result.

| Region | Primary 850K rho | Coordinate-WGBS rho | WGBS neuron FDR | WGBS pass |
|---|---:|---:|---:|:---:|
| DLPFC | -0.035 | -0.020 | 0.798 | No |
| Hippocampus | -0.103 | -0.107 | 0.317 | No |

Lee- and Tian-specific coordinate-WGBS estimates also failed the neuronal gate.
Thus, the alternative reference does not rescue either region and strengthens
the decision not to use their DNAm fractions as meQTL covariates.

## Decision

- Retain all three regions in the DNAm deconvolution QC report and figure.
- Permit M6d only for caudate.
- Build `M3a_dnam_matched` and M6d on the same 151 AA caudate samples and the
  same CpG phenotype BED; do not impute the two failed AA samples.
- Do not build M6d for DLPFC or hippocampus.
- Do not promote M6d into the primary analysis. Because fewer than two regions
  pass integration QC, the prespecified multi-region manuscript-retention
  criterion cannot be satisfied even before TensorQTL sensitivity mapping.
- RNA MuSiC remains an orthogonal validation source and its files are unchanged.

The final gate values are machine-readable in
`_m/dnam-scmd-validation-summary.tsv`.
