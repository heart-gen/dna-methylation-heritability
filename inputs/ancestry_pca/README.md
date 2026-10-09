# inputs/ancestry_pca

Places the cohort in the 1000 Genomes principal-component space, for the
ancestry-confirmation supplement (`figureS_ancestry_pcs`, Module 11).

One PCA is fit on the 3,202 1000 Genomes samples alone, with allele weights.
Both the reference samples and the cohort's donors are then scored with those
same weights, so every point shares one basis. An earlier version overlaid the
cohort's own PCs on a separately computed 1000 Genomes PCA; two PCAs have
unrelated axes, and the overlay put white American donors beside AFR.

This is a QC panel. It confirms the recorded donor group against genotype and
supports no ancestry-biology claim.

## Variants

Autosomal, biallelic A/C/G/T SNPs with MAF >= 0.05 and missingness <= 0.05 in
both panels, keyed on chr:pos, with identical allele sets and no
strand-ambiguous (A/T, C/G) pairs, LD-pruned in the reference
(`--indep-pairwise 500 50 0.1`). `_m/variant_counts.tsv` records each step.

## Run

```bash
cd inputs/ancestry_pca/_m && mkdir -p logs
sbatch ../_h/step_1_project_1kgp.sh
```

Reads `inputs/genotypes/all_individuals/TOPMed_LIBD` and the 1000 Genomes
pfile under `/projects/b1213/resources/1kGP/GRCh38_phased_vcf/_m/`.

## Outputs (`_m/`, tracked)

| file | content |
|---|---|
| `projected_pcs.tsv` | `source` (1000 Genomes / this study), `id` (1000 Genomes sample or BrNum), PC1-PC10 |
| `1kgp_shared.eigenval` | eigenvalues of the reference PCA |
| `variant_counts.tsv` | variants at each filtering step |
| `output_checksums.tsv` | SHA-256 of the three files above |

Donor identifiers are BrNum; the array barcode is dropped.
