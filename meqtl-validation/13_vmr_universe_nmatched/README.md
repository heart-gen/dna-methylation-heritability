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
step 1 creates essentially the entire difference.

The "mean sd cutoff" column is the pipeline's own stored value and is shown for
completeness only. **Do not read a between-region variability claim off it** —
it does not reproduce, and the cutoffs are not comparable across regions because
DLPFC and hippocampus computed theirs on a truncated donor set. See
[Open item](#open-item-reconcile-pccsv-before-any-cross-region-sd-claim).

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

Full array completed 2026-08-15 (66/66 tasks, 0 errors). **The caudate universe
surplus is neither sample size nor biology — it is library coverage uniformity.**

### Genome-wide CpG universe (autosomes; same 25,802,333 candidate sites per region)

| Design | N | CpGs passing QC | |
|---|---:|---:|---|
| Caudate, full N | 153 | 23,466,905 (90.9%) | |
| **Caudate, N-matched** | **111** | **23,476,406** | mean of 30 reps; range 23,431,833–23,530,379 |
| DLPFC | 111 | 20,521,441 (79.5%) | |
| Hippocampus | 116 | 20,764,810 (80.5%) | |

### Comparison ladder

| Contrast | Controls for | Ratio |
|---|---|---:|
| caudate N-matched vs caudate full-N | everything but donor count | **1.0004** |
| caudate N-matched vs DLPFC | donor count | 1.1440 |
| caudate N-matched vs hippocampus | donor count | 1.1306 |

**Sample size explains 0.04% of the surplus.** Dropping 42 donors moves the
caudate universe by ~9,500 CpGs out of 23.5M, and the full-N value sits inside
the replicate range. The ~14% surplus over both comparators survives intact.

### Why, if not sample size — and it is not depth either

| Region | Mean cov | SD cov | Frac sites cov≥5 | SD frac | Min frac |
|---|---:|---:|---:|---:|---:|
| Caudate | 17.47 | **1.93** | **0.9216** | **0.0307** | 0.680 |
| DLPFC | 18.27 | 5.08 | 0.8582 | 0.0640 | 0.492 |
| Hippocampus | 18.24 | 5.55 | 0.8670 | 0.0557 | 0.698 |

Paired on the 92 donors assayed in all three regions:

| Contrast | Δ depth | Δ adequacy |
|---|---:|---:|
| caudate − DLPFC | **−0.665** | +0.0641 |
| caudate − hippocampus | **−0.710** | +0.0544 |

Caudate is **shallower** than both comparators on the same donors, yet keeps
more CpGs. The resolution is dispersion: caudate's between-donor coverage SD is
1.93 against 5.08 and 5.55, and its per-donor adequate-site fraction is both
higher and ~2× tighter.

This follows from the filter's form. `exclude_low_cov()` keeps a CpG when **≥80%
of donors** clear cov≥5 — a property of the low tail of the donor distribution,
not of mean depth. A region with a few poorly covered donors (DLPFC has one at
0.49 adequacy) loses CpGs genome-wide however deep its good libraries are.

### The variability bar: within-caudate robust, between-region not

Caudate's own bar barely moves under N-matching — 0.070605 full-N vs 0.070144
N-matched. Same method on both sides, so this contrast is trustworthy: the
relative-quantile cutoff is **not** an N artifact.

The between-region comparison is **not robust and should not be used**:

| Region | Module 13 (recomputed) | Pipeline `top1_cpg.tsv` |
|---|---:|---:|
| Caudate | 0.070605 | 0.079675 |
| DLPFC | 0.065246 | 0.080218 |
| Hippocampus | 0.067412 | 0.081225 |

The two disagree on the ordering — caudate has the highest bar on 17/22
autosomes here but only 7/22 in the pipeline's own output, and the two series
correlate at only r=0.377 across 66 region×chromosome cells. This module
reproduces the pipeline's two-stage adjustment (regress on snpPC1–3, PCA the
residuals, regress on 5 methylation PCs), but estimates the PCs differently
from `02.pca.R`. So **no claim about which region's VMRs sit at a higher
variability bar is supportable** from either source until that is reconciled.

Tracing the discrepancy surfaced a defect in `02.pca.R` that makes the
pipeline's own column the less trustworthy of the two. See below.

## Open item: reconcile `pc.csv` before any cross-region sd claim

`02.pca.R` writes the methylation PCs (`pc.csv`, per region × chromosome) that
`02b.res_var.R` regresses out before the top-1% VMR seed cutoff. Four
differences separate it from this module, one of which is a bug.

**1. `get_snp_pcs()` hardcodes `region == "caudate"` in all three regions'
copies.** `vmr-analysis/{caudate,dlpfc,hippocampus}/_h/02.pca.R:34` is identical
in each. snpPCs are per-donor so the *values* are right, but the filter
restricts the donor set to donors having a caudate sample, and `res_snp_pcs()`
then subsets methylation to those survivors:

| Region | Donors w/ snpPCs | Also have caudate | **Dropped** |
|---|---:|---:|---:|
| Caudate | 153 | 153 | **0** |
| DLPFC | 111 | 96 | **15** |
| Hippocampus | 116 | 101 | **15** |

`02b.res_var.R` then does `valid_ids <- intersect(valid_ids, pc$V1)`, so the
loss propagates: DLPFC and hippocampus residual SDs — and therefore their
top-1% cutoffs, VMR seeds, and VMR intervals — were computed on **96 and 101
donors**, while caudate used all 153. All 30 dropped donors reappear in the
downstream `_m/samples.txt` (111/116), so the regions' VMRs were *defined* on a
truncated donor set but *modeled* on the full one.

The asymmetry runs in the direction under dispute: it is caudate-sparing, and
fewer donors gives a noisier SD estimate and a higher 99th percentile — which is
what the pipeline column shows for DLPFC and hippocampus.

**2. CpG selection.** `get_top_meth()` takes the top 1M CpGs by *raw*
(unresidualized) SD; this module uses a random 50,000 coverage-passing CpGs.
Variance-enriched and deterministic versus unbiased and random.

**3. Scaling.** The pipeline uses `prcomp(scale. = TRUE)`, this module
`scale. = FALSE`. On variance-ranked input, unit-scaling changes what PC1–5
absorb.

**4. Site count.** Up to 10^6 versus 5×10^4.

Also cosmetic, in the DLPFC and hippocampus copies only: `v_top[1:10^6, ]`
without caudate's `min(10^6, nrow(v_top))` guard, which pads with NA rows on
chromosomes holding under 1M CpGs. Harmless — the NAs never match in
`is.element()` — but it should carry the same guard.

### To do

1. Fix the hardcoded region and re-run `02.pca.R` → `02b` → `02c` for DLPFC and
   hippocampus at their true N. This is required regardless of the
   reconciliation; it changes the published VMR sets.
2. Re-diff the two sd series after the fix, then vary selection, scaling, and
   site count one at a time to close whatever gap remains.
3. Only then decide whether a cross-region variability statement is supportable.

Nothing in this module's headline result depends on the outcome. The
`coverage_uniformity_not_sample_size` verdict rests on CpG counts under the
coverage filter, which never touches `pc.csv`, and the within-caudate contrast
(0.070605 → 0.070144) uses one method on both sides in a region with zero
dropped donors.

**Not yet verified empirically:** `heritability/<region>/_m/pca/` does not exist
in this repo, so no `pc.csv` survives to confirm the row counts. The defect is
read from the code and the phenotype table. Alexis's repo has caudate-only
copies under `covar-analysis/caudate/{snp_pc,no_snp_pc,global_ances}/_m/pca/`.

## What this means for the caudate claim

"Caudate's greater discovery is not solely sample size" is **literally
supported** — sample size accounts for 0.04%. But the mechanism that does
account for it is technical, not biological: caudate contributes ~14% more
testable CpGs because its WGBS libraries are more uniform across donors.

So module 10's common-universe restriction is doing the right thing. It removes
CpGs testable only in caudate, and those are caudate-only for a library-
uniformity reason. The failed gate should be read neither as "caudate is
unremarkable" nor as "the restriction is over-conservative" — the honest
statement is that the raw-count advantage is a technical artifact and the
common-universe rates (0.4337 vs 0.4601 / 0.4674) are the comparable numbers.

Not addressed here: λ_GC = 3.884 in module 10.

## Data-integrity issues found and resolved

**Broken HDF5 backing (resolved).** The DLPFC and hippocampus
`<region>_chr<N>_BSobj.rda` archives in this repo cannot be loaded: `M`/`Cov`
point at `combined_hdf5/<region>_assays.h5` where the file on disk is
`combined_hdf5/assays.h5`, and `coef` points at a deleted R session temp dump.
Alexis re-exported both regions with `saveHDF5SummarizedExperiment()` — a
self-contained `assays.h5` + `se.rds` directory with all three assays resolving
and `hasBeenSmoothed()` TRUE, covering all 22 autosomes. `load_bsobj()` prefers
those and falls back to the local `.rda` for caudate, whose backing is intact.
Coverage results are identical either way, confirming the re-export is the same
data. Without it the SD arm could only have run on caudate.

**Donor count off-by-one (resolved).** Earlier runs found 154/112/117 donors
against module 10's design of 153/111/116. The extra donor is `Br1442`, the same
one in all three regions, who has no genotype PCs. `02.pca.R::get_snp_pcs()`
drops NA-PC samples and `res_snp_pcs()` subsets methylation to the survivors, so
the pipeline never included them. This module now applies the same filter. It is
not just bookkeeping — the coverage filter asks whether ≥80% of *these* donors
clear cov≥5, so an ungenotyped donor shifts the universe.

## Run

```bash
cd meqtl-validation/13_vmr_universe_nmatched/_m
mkdir -p logs
sbatch ../_h/step_1.sh    # 66 tasks: 22 autosomes x 3 regions, ~4 min each
/projects/p32505/opt/envs/genomics/bin/python3 ../_h/02_aggregate_universe.py --require-complete
```

## Outputs

Per region × chromosome (66 each):

- `_m/coverage_pass.{region}.chr{N}.tsv` — CpGs passing the filter, full-N and per replicate
- `_m/donor_coverage.{region}.chr{N}.tsv` — per-donor mean coverage and adequate-site fraction
- `_m/sd_cutoff.{region}.chr{N}.tsv` — residual-sd 99th percentile

Aggregated:

- `_m/universe_totals.tsv` — genome-wide pass counts per region and replicate
- `_m/universe_comparisons.tsv` — the comparison ladder
- `_m/sd_cutoff_summary.tsv` — variability bar by region and design
- `_m/universe_verdict.tsv` — machine-readable call, `coverage_uniformity_not_sample_size`
