# Module 13: Is the caudate VMR/CpG universe advantage sample size?

## Question

Module 10 downsampled donors for meQTL **mapping** but inherited the VMR/CpG
universe from full-N caudate. The universe is therefore still confounded — and
it is the universe, not the mapping, that drives the common-universe rate
comparison that failed the caudate gate.

Caudate contributes 11,373 VMRs / 197,466 CpGs against DLPFC 9,976 / 154,325
and hippocampus 9,801 / 156,665, at N=153 vs 111 / 116.

## Where the difference is created

Two steps in `vmr-analysis/<region>/_h` build the universe:

1. `01.get_cpg_stats.R::exclude_low_cov()` keeps CpGs with
   `rowSums2(coverage >= 5) >= n_donors * 0.8`.
2. `02c.write_top_cpg.R` takes VMR seeds as CpGs above the **within-region,
   per-chromosome 99th percentile** of residual sd — a *relative* cutoff, so
   seed count is mechanically ~1% of whatever survives step 1.

Measured across autosomes at full N:

| Region | CpGs passing QC | Top-1% seeds | Mean sd cutoff |
|---|---:|---:|---:|
| Caudate | 24,435,609 | 244,367 | 0.0822 |
| DLPFC | 21,115,720 | 211,168 | 0.0854 |
| Hippocampus | 21,383,452 | 213,845 | 0.0865 |

Caudate keeps **15.7% more CpGs** and therefore gets 15.7% more seeds by
construction, which tracks the VMR surplus (+14% / +16%) almost exactly. So
step 1 creates essentially the entire difference, and caudate's "top 1%" sits
at a **lower absolute variability bar** than the comparators'.

## Design

Re-run only the coverage filter on the same 30 N=111 caudate donor subsets
module 10 used (seed 20260805), against full-N caudate and both comparators.
Also measure per-donor coverage directly, and recompute the residual-sd 99th
percentile on each universe.

- caudate N-matched ≈ caudate full-N, still above comparators → **coverage depth**
- caudate N-matched falls to comparator level → **sample-size artifact**

## Run

```bash
cd meqtl-validation/13_vmr_universe_nmatched/_m
mkdir -p logs
sbatch ../_h/step_1.sh    # 66 tasks: 22 autosomes x 3 regions
```

## Result

Full array completed 2026-08-15 (66/66 tasks). **The caudate universe surplus is
neither sample size nor biology — it is library coverage uniformity.**

### Genome-wide CpG universe (autosomes; same 25,802,333 candidate sites in every region)

| Design | N | CpGs passing QC | vs comparator |
|---|---:|---:|---|
| Caudate, full N | 153 | 23,451,375 (90.9%) | — |
| **Caudate, N-matched** | **111** | **23,476,406** | mean of 30 reps; range 23,431,833–23,530,379 |
| DLPFC | 111 | 20,507,140 (79.5%) | — |
| Hippocampus | 116 | 20,757,353 (80.4%) | — |

### Comparison ladder

| Contrast | Controls for | Ratio |
|---|---|---:|
| caudate N-matched vs caudate full-N | everything but donor count | **1.0011** |
| caudate N-matched vs DLPFC | donor count | 1.1448 |
| caudate N-matched vs hippocampus | donor count | 1.1310 |

**Sample size explains 0.1% of the surplus.** Dropping 42 donors moves the
caudate universe by +25,031 CpGs out of 23.5M, and the full-N value sits inside
the replicate range. The ~14% surplus over both comparators survives N-matching
intact.

### Why, if not sample size — and it is not depth either

| Region | Mean cov | SD cov | Frac sites cov≥5 | SD frac | Min frac |
|---|---:|---:|---:|---:|---:|
| Caudate | 17.45 | **1.95** | **0.9214** | **0.0307** | 0.6801 |
| DLPFC | 18.25 | 5.06 | 0.8582 | 0.0637 | 0.4921 |
| Hippocampus | 18.23 | 5.53 | 0.8672 | 0.0555 | 0.6984 |

Paired on the 92 donors assayed in all three regions:

| Contrast | Δ depth | Δ adequacy |
|---|---:|---:|
| caudate − DLPFC | **−0.665** | +0.0641 |
| caudate − hippocampus | **−0.710** | +0.0544 |

Caudate is **shallower** than both comparators on the same donors, yet keeps
more CpGs. The resolution is dispersion: caudate's between-donor SD of coverage
is 1.95 against 5.06 and 5.53, and its per-donor adequate-site fraction is both
higher and ~2× tighter.

This follows directly from the filter's form. `exclude_low_cov()` keeps a CpG
when **≥80% of donors** clear cov≥5 — a statement about the low tail of the
donor distribution, not about mean depth. A region with a handful of poorly
covered donors (DLPFC has one at 0.49) loses CpGs genome-wide regardless of how
deep its good libraries are. Caudate's libraries are more consistent, so it
loses fewer.

### The variability bar barely moves under N-matching

Caudate's own residual-sd 99th percentile: 0.071563 at full N, 0.071413
N-matched (22 autosomes, 5 replicates each). So the relative-quantile cutoff is
not itself an N artifact.

## What this means for the caudate claim

The claim "caudate's greater discovery is not solely sample size" is **literally
supported** — sample size accounts for 0.1%. But the mechanism that does account
for it is technical, not biological: caudate contributes ~14% more testable CpGs
because its WGBS libraries are more uniform across donors.

So module 10's common-universe restriction is doing the right thing. It removes
CpGs that are testable only in caudate, and those CpGs are caudate-only for a
library-uniformity reason. The failed gate should not be read as "caudate is
unremarkable" nor defended as "the restriction is over-conservative" — the
honest statement is that the raw-count advantage is a technical artifact, and
the common-universe rates (0.4337 vs 0.4601 / 0.4674) are the comparable numbers.

Not addressed here: λ_GC = 3.884 in module 10.

## Caveats

- The sd-cutoff arm is **caudate only**. It needs smoothed methylation, and the
  DLPFC/hippocampus smoothing coefficients are unrecoverable (see below). The
  cross-region variability bars quoted in the section above come from the
  pipeline's own `vmr-analysis/<region>/_m/cpg/top1_cpg.tsv`.
- The DLPFC and hippocampus `BSobj` archives have broken HDF5 backing. `M`/`Cov`
  point at `<region>_assays.h5` where the file on disk is `assays.h5` — a
  recoverable rename, repointed at load time. `coef` points at a deleted R
  session temp dump (`auto<hash>.h5`) and is gone. `01_universe_decomposition.R`
  repairs the former and drops the latter; it does not modify the shared
  archives, which are owned by another user.
- Donor counts found here (154 / 112 / 117) exceed module 10's design summary
  (153 / 111 / 116) by one in each region. Immaterial to these ratios, but worth
  resolving before the numbers reach a figure.

## Status

Complete. Aggregate outputs regenerate in seconds:

```bash
cd meqtl-validation/13_vmr_universe_nmatched/_m
/projects/p32505/opt/envs/genomics/bin/python3 ../_h/02_aggregate_universe.py
```

## Outputs

Per region x chromosome (66 each):

- `_m/coverage_pass.{region}.chr{N}.tsv` — CpGs passing the filter, full-N and per replicate
- `_m/donor_coverage.{region}.chr{N}.tsv` — per-donor mean coverage and adequate-site fraction
- `_m/sd_cutoff.{region}.chr{N}.tsv` — residual-sd 99th percentile (caudate only)

Aggregated:

- `_m/universe_totals.tsv` — genome-wide pass counts per region and replicate
- `_m/universe_comparisons.tsv` — the comparison ladder
- `_m/sd_cutoff_summary.tsv` — variability bar by design
- `_m/universe_verdict.tsv` — machine-readable call, `coverage_uniformity_not_sample_size`
