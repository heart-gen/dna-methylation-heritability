# Phase: official caudate TensorQTL downsampling

**Status: Pass / complete** (30/30 replicates; claim snapshot written).

N-matched caudate cis-meQTL remapping with **TensorQTL permutation / beta-approximated FDR** (same engine as Phase 1), using the Experiment 3 shared-donor-aware sample lists (target N=111, 30 replicates).

This supersedes lead-SNP retention absolute discovery-count comparisons for the formal caudate≠N claim. Lead-SNP retention (`04_cross_region_sharing/_m/caudate_downsample/`) remains a fast sensitivity.

## Result snapshot

| Metric | Value |
|---|---:|
| Full caudate TensorQTL FDR sig CpGs | 81,007 |
| Downsample median n_sig (N=111) | **67,562.5** (IQR 67,181–67,922) |
| Median fraction of full-N FDR hits retained | **0.783** |
| DLPFC / hippocampus Phase 1 n_sig | 51,273 / 55,957 |
| Median ratio vs DLPFC / hip | 1.32 / 1.21 |
| `criterion_not_solely_sample_size` | **True** (all 30 reps exceed both comparators) |

Source: `_m/tensorqtl_downsample_claim_snapshot.tsv`

## Design

| Item | Value |
|---|---|
| Region | caudate (AA) |
| Covariates | Locked **M3a** (`01_.../caudate/_m/prepared/covariates.txt`) |
| Genotypes | `01_.../caudate/_m/genotypes/meqtl_AA` |
| Phenotypes | Subset of `prepared/cpg_phenotypes.all_autosomes.bed.gz` |
| Sample lists | Reused from `04_cross_region_sharing/_m/caudate_downsample/sample_lists/` (seed 20260805) |
| Target N | 111 (= min DLPFC, hippocampus) |
| Replicates | 30 |
| Cis window / MAF / FDR | ±500 kb / 0.05 / 0.05 |

## Run

```bash
cd meqtl-validation/10_downsampling_caudate/_m
mkdir -p logs

# CPU: copy lists + write per-rep covariates + phenotype BEDs
J1=$(sbatch --parsable ../_h/step_1.sh)

# GPU array: official TensorQTL cis + QC per replicate
J2=$(sbatch --parsable --dependency=afterok:${J1} ../_h/step_2.sh)

# CPU: aggregate vs full caudate / DLPFC / hippocampus
sbatch --dependency=afterok:${J2} ../_h/step_3.sh
```

Optional smoke test (1 rep):

```bash
sbatch --export=ALL,MAX_REPS=1 ../_h/step_1.sh
sbatch --array=0 --export=ALL,MAX_REPS=1 ../_h/step_2.sh
```

## Outputs

```text
_m/
  sample_lists/                 # symlinks or copies of Exp3 lists
  downsample_design_summary.tsv
  downsample_replicate_manifest.tsv
  rep001/ ... rep030/
    covariates.txt
    cpg_phenotypes.bed.gz (+ .tbi)
    tensorqtl/*.cis_qtl.txt.gz
    tensorqtl/qc/meqtl_qc_summary.tsv
  tensorqtl_downsample_replicate_results.tsv
  tensorqtl_downsample_vs_regions.tsv
  tensorqtl_downsample_claim_snapshot.tsv
```

## Interpretation rule

Compare **TensorQTL permutation-FDR** n_sig from downsampled caudate to DLPFC/hippocampus Phase 1 M3a n_sig (same method family). Prefer median n_sig / retention of full-N FDR hits over lead-SNP BH-FDR counts from module 04.
