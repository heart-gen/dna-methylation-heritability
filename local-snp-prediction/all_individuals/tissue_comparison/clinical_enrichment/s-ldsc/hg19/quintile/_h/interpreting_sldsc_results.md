# Interpreting S-LDSC Results

## Overview

Stratified LD Score Regression (S-LDSC) partitions SNP-heritability across
genomic annotations to test whether specific categories of SNPs are enriched
for trait heritability. Unlike basic LDSC, which outputs a single genome-wide
h2 estimate, S-LDSC produces per-annotation enrichment statistics.

## Output Files

Each S-LDSC run produces two key files:

### `.log` file

Contains the **total observed-scale SNP-heritability** estimate, e.g.:

```
Total Observed scale h2: 0.0761 (0.0166)
```

This is the single h2 value (with standard error) for the GWAS trait. It also
reports:

| Field | Description |
|---|---|
| **Lambda GC** | Genomic inflation factor; values close to 1 suggest minimal confounding |
| **Mean Chi^2** | Average chi-squared statistic; should be > 1 if there is signal |
| **Intercept** | LDSC intercept (SE); values near 1 indicate polygenicity rather than confounding |
| **Ratio** | (Intercept - 1) / (Mean Chi^2 - 1); proportion of inflation due to confounding vs. polygenicity |

### `.results` file

A tab-separated table with one row per annotation. Each row reports how much
of the trait's heritability is attributable to SNPs in that annotation.

## Column Definitions

| Column | Definition | Interpretation |
|---|---|---|
| **Category** | Annotation name | Genomic feature (e.g., enhancer, promoter, coding region) |
| **Prop._SNPs** | Proportion of reference SNPs in this annotation | Size of the annotation |
| **Prop._h2** | Proportion of h2 explained by this annotation | Fraction of heritability captured |
| **Prop._h2_std_error** | Standard error of Prop._h2 | Uncertainty in the h2 proportion |
| **Enrichment** | Prop._h2 / Prop._SNPs | Fold-enrichment; values > 1 indicate the annotation captures more h2 than expected by its size |
| **Enrichment_std_error** | Standard error of enrichment | Uncertainty in the enrichment estimate |
| **Enrichment_p** | P-value for enrichment | Whether the enrichment is statistically significant |
| **Coefficient** | Regression coefficient (tau) | Per-SNP contribution to h2 **after controlling for all other annotations** |
| **Coefficient_std_error** | Standard error of tau | Uncertainty in the coefficient |
| **Coefficient_z-score** | tau / SE(tau) | Significance of the annotation's unique contribution; |z| > 1.96 corresponds to p < 0.05 |

## Key Metrics to Report

For most analyses, focus on:

1. **Enrichment + Enrichment_p** — Is the annotation enriched for heritability?
   A significant enrichment (p < 0.05) with fold-enrichment > 1 means SNPs in
   this annotation explain more heritability than expected given the number of
   SNPs.

2. **Coefficient_z-score** — Does the annotation contribute to heritability
   *beyond* what is explained by overlapping baseline annotations? This is the
   more stringent test because it conditions on all other annotations in the
   model.

## Identifying the Custom Annotation

The `.results` file contains rows for all annotations in the model:

- **Rows labeled `*L2_0`** — These are the 97 baseline LD model annotations
  (v2.2) from the Broad Institute, including functional categories (coding,
  enhancer, promoter, etc.), histone marks, conservation scores, MAF bins, and
  other genomic features.

- **The last row (`L2_1`)** — This is the **custom annotation** added on top
  of the baseline model. In this analysis, `L2_1` corresponds to the
  tissue-specific CpG annotation (e.g., caudate heritable CpGs). This is the
  row of primary interest for testing whether DNA methylation heritability
  categories are enriched for GWAS trait heritability.

### Example interpretation of L2_1

```
Category    Prop._SNPs  Prop._h2   Enrichment  Enrichment_p  Coefficient_z-score
L2_1        0.0043      -0.0305    -7.04       0.178         -1.37
```

- The custom annotation contains ~0.43% of SNPs
- Enrichment is negative (no enrichment for trait heritability)
- Enrichment p-value is non-significant (p = 0.178)
- Coefficient z-score is non-significant (|z| < 1.96)
- **Conclusion**: This annotation is not enriched for the trait's heritability

## Extracting Results Across Traits

To extract the total h2 from all log files:

```bash
grep "Total Observed scale h2" results/*/*/*/*.log
```

To extract only the custom annotation row (L2_1) from all results files:

```bash
# Header
head -1 results/ad/caudate/heritable_hg19/ad_caudate_heritable_hg19.results

# Custom annotation row from all results
tail -1 results/*/*/*/*.results
```

## References

- Bulik-Sullivan et al. (2015). LD Score regression distinguishes confounding
  from polygenicity in genome-wide association studies. *Nature Genetics*.
- Finucane et al. (2015). Partitioning heritability by functional annotation
  using genome-wide association summary statistics. *Nature Genetics*.
- Gazal et al. (2017). Linkage disequilibrium-dependent architecture of human
  complex traits shows action of negative selection. *Nature Genetics*.
