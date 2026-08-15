# VMR pipeline audit — defect register and task list

Audit date 2026-08-15. Scope: `vmr-analysis/**/_h/*` (both the `BA_only` branch at
`vmr-analysis/<region>/` and the `all_individuals` branch), checked against the
production outputs in the lead author's repo,
`/projects/b1213/users/alexis/projects/dna-methylation-heritability`.

**All six VMR scripts are byte-identical between the two repos** for all three
regions, and `_m/cpg/top1_cpg.tsv` is identical as well. The two repos are the
same analysis; Alexis's copy additionally retains the intermediates
(`_m/pca/`, `_m/cpg/chr_*/tmp_files/`, `_m/vmr/`) that this repo does not, which
is what made empirical verification possible.

Seven defects, one of which invalidates every VMR set in the project.

---

## V1 — CRITICAL: `02b.res_var.R` regresses methylation on permuted donor rows

`vmr-analysis/*/_h/02b.res_var.R` and `vmr-analysis/all_individuals/*/_h/02b.res_var.R`,
in `filter_pheno()` (BA_only lines 26–27):

```r
meth_levels <- meth_levels[match(valid_ids, brain_id), , drop = FALSE]  # valid_ids order
pc_filt     <- pc[V1 %in% valid_ids, -1, with = FALSE]                  # pc.csv row order
```

`valid_ids` comes from `intersect(intersect(ances$id, demo$brnum), brain_id)`, so it
carries the **ancestry file's** row order. `pc[V1 %in% valid_ids]` is a data.table
logical subset, so it preserves **`pc.csv`'s** row order. The two orders are unrelated.
`regress_pcs_vectorized()` then builds `design <- cbind(1, pc_matrix[, 1:5])` and runs
`lmFit(t(meth_levels), design)` — response rows and design rows belong to different
donors.

Measured misalignment at chr1: caudate **153/153** rows wrong, DLPFC 92/96,
hippocampus 98/101.

### Verification

Recomputing residual SDs for caudate chr1 chunk `cpg_meth_10003_15002.tsv` both ways
and comparing against the pipeline's stored `_m/pca/chr_1/res_var_all.tsv`:

| Version | Pearson r vs stored | max abs diff |
|---|---:|---:|
| Correctly aligned | 0.991245 | 3.1×10⁻² |
| **As coded (pc.csv order)** | **1.000000** | **5.0×10⁻¹⁶** |

The as-coded version reproduces the shipped output to machine precision. This is not
a hypothesis.

### Impact

Full recomputation of caudate chr21 (336,260 CpGs, 153 donors), correct vs as-coded:

| Quantity | Correct | As coded | |
|---|---:|---:|---|
| 99th-pct sd cutoff | 0.081376 | 0.088864 | **+9.2%** |
| Top-1% seed CpGs | 3,363 | 3,363 | fixed by construction |
| Seed overlap | — | 3,032 | **Jaccard 0.82** |

About **10% of VMR seed CpGs are wrong**, and the variability bar is inflated ~9%
because permuted PCs remove noise instead of real structure. Since VMRs require ≥6
consecutive high CpGs within 1 kb, VMR-level turnover is not the same as seed
turnover and has not been measured.

Affects **all 6 script copies** — both branches, all three regions — so every VMR set
in the project is derived from a mis-specified regression, including the `vmr.bed`
files consumed by the elastic-net heritability analysis and all of `meqtl-validation`.

**Fix:** order the PC matrix by `valid_ids`, e.g.
`pc_filt <- pc[match(valid_ids, V1), -1, with = FALSE]`, and add a hard
`stopifnot(identical(pc_ids, valid_ids))` guard.

---

## V2 — HIGH: `02.pca.R` hardcodes `region == "caudate"` in all three BA_only copies

`vmr-analysis/{caudate,dlpfc,hippocampus}/_h/02.pca.R:34` is identical in all three:

```r
snp_pcs <- fread(pc_file_path, header=TRUE)[
    race == "AA" & agedeath >= 17 & region == "caudate", ...]
```

snpPCs are per-donor so the values are right, but the filter restricts the donor set
to donors who have a caudate sample, and `res_snp_pcs()` subsets methylation to the
survivors.

**Verified in the production outputs** — `_m/pca/chr_*/pc.csv` row counts, uniform
across all 22 autosomes:

| Region | Expected | Actual | Dropped |
|---|---:|---:|---:|
| Caudate | 153 | 153 | 0 |
| DLPFC | 111 | **96** | **15** |
| Hippocampus | 116 | **101** | **15** |

`02b.res_var.R` then does `valid_ids <- intersect(valid_ids, pc$V1)`, so the loss
propagates into the residual SDs, the top-1% cutoff, and the VMR intervals. But
`04.cal_vmr.R` loads `stats.rda`, whose `BSobj` still holds 111/116 donors — so
DLPFC and hippocampus VMRs were **defined** on 96/101 donors and **modeled** on
111/116. `_m/samples.txt` confirms 153/111/116 downstream.

The `all_individuals` copies have the correct per-region filter (`pc.csv` = 282/166/176),
which confirms this is a copy-paste defect rather than an intended restriction. It is
also caudate-sparing, so it biases the cross-region comparison in the direction of the
disputed caudate gate, and it is the leading explanation for module 13's failure to
reproduce the between-region variability bar (r = 0.377).

---

## V3 — HIGH: committed scripts do not reproduce committed results (path drift)

Commit `5a6ae1570` (2026-02-05, "Move vmr pipeline to separate dir") renamed
`heritability/<region>/` → `vmr-analysis/<region>/` as a pure R100 rename. The
`here()` paths *inside* the scripts were only partly updated.

| Script | Reads/writes | Directory exists? |
|---|---|---|
| `01.get_cpg_stats.R` (all regions) | `heritability/<region>/_m` | **no** |
| `02.pca.R` (all regions) | `heritability/<region>/_m/pca` | **no** |
| `02b.res_var.R` (all regions) | `heritability/<region>/_m/{cpg,pca}` | **no** |
| `02c.write_top_cpg.R` (all regions) | `vmr-analysis/<region>/_m/{pca,cpg}` | yes |
| `03.extract_vmr.R` (caudate) | `vmr-analysis/caudate/_m/{pca,vmr}` | yes |
| `03.extract_vmr.R` (DLPFC, hippocampus) | `heritability/<region>/_m/{pca,vmr}` | **no** |
| `04.cal_vmr.R` (all regions) | `heritability/<region>/_m/{cpg,vmr}` | **no** |

`heritability/` now contains only `elastic_net_model/` in both repos. The pipeline
cannot run end-to-end as committed, and the stage boundary between `02b` (writes
`heritability/`) and `02c`/`03` (reads `vmr-analysis/`) is broken. This must be fixed
before any re-run, so it gates V1 and V2.

---

## V4 — MEDIUM: sex chromosomes are in the VMR sets and are not C→T masked

`step_3.sh` is `--array=1-24`, and `02c.write_top_cpg.R` iterates `c(1:22, "X", "Y")`,
so VMRs are called on X and Y. But `01.get_cpg_stats.R:134` skips CT-SNP removal for
them:

```r
if (!chr %in% c("X", "Y")) { BSobj <- remove_ct_snps(f_snp, filtered$BSobj) }
```

because `/projects/b1213/resources/libd_data/wgbs/DEM2/snps_CT/` only holds `chr1`–`chr22`.
So sex-chromosome VMRs are the only ones contaminated by C→T SNPs, and they sit in a
region where X-inactivation and sex dosage inflate between-donor variance.

| Region | Total VMRs | X/Y VMRs | % female |
|---|---:|---:|---:|
| Caudate | 12,001 | **431** | 35% |
| DLPFC | 10,372 | 143 | 41% |
| Hippocampus | 10,216 | 147 | 41% |

Caudate carries a **3× excess** of sex-chromosome VMRs that sex composition does not
explain and that is far larger than its 15.7% autosomal universe surplus. The cause is
not diagnosed. These are not inert: `heritability/elastic_net_model/*/01.elastic-net.R`
recodes X→23 and Y→24 and carries them into the predictability analysis.

Decide whether to drop X/Y from VMR calling, or to mask and justify them — and
diagnose the caudate excess either way.

---

## V5 — LOW: `remove_ct_snps()` ignores its own argument

`01.get_cpg_stats.R:30` reads the global `filtered` rather than the `BSobj` parameter:

```r
remove_ct_snps <- function(f_snp, BSobj) {
    idx <- is.element(start(filtered$BSobj), snp)   # global, not the argument
    BSobj <- BSobj[!idx, ]
```

Harmless today because the only call site is `remove_ct_snps(f_snp, filtered$BSobj)`,
so the two coincide. Silently wrong for any other call, and it defeats the guard the
parameter was meant to provide.

## V6 — LOW: missing `nrow` guard in the DLPFC and hippocampus `02.pca.R`

`02.pca.R:20` is `v_top[1:10^6, ]` in those two copies but `v_top[1:min(10^6, nrow(v_top)), ]`
in caudate's. On chromosomes holding under 1M CpGs the un-guarded form pads `v_top`
with NA rows. Benign — the NAs never match in `is.element()` — but the guard should be
uniform.

## V7 — LOW: `res_cpg_meth.phen` columns are in lexicographic, not genomic, order

`02b.res_var.R:134` collects chunks with `list.files()`, which sorts lexicographically
(`residuals_1000003_…` before `residuals_100003_…`), then `cbind`s them. Columns are
named so the data is recoverable, but any positional read is wrong. The `cbind` also
assumes every chunk shares the first chunk's FID/IID order without checking. Only
consumer today is `meqtl-validation/00_data_audit/_h/00_data_audit.py`.

---

## Checked and cleared

- `quantile(v[, 3], ...)` in `03.extract_vmr.R` vs `quantile(v$sd, ...)` in
  `02c.write_top_cpg.R` — tested in R; the one-column data.table gives the identical
  result. Not a defect.
- Region filters in `01.get_cpg_stats.R` and `02b.res_var.R` are correct per region.
  Only `02.pca.R` carries the hardcode.
- All six VMR scripts and `top1_cpg.tsv` are identical between this repo and Alexis's.

---

## Task list

Ordering is forced by data dependency: V3 unblocks any re-run, V1 and V2 must land in
the same re-run, and everything downstream of VMR definition has to be re-verified
afterward.

### Stage A — repair and re-run the VMR pipeline

| # | Task | Depends on |
|---|---|---|
| **V3** | Repoint every `here()` in `01`, `02`, `02b`, `04` (all regions) and `03` (DLPFC, hippocampus) from `heritability/<region>` to `vmr-analysis/<region>`. Add a smoke test that the stage boundaries resolve. | — |
| **V1** | Fix the PC row misalignment in all 6 `02b.res_var.R` copies; add an `identical()` guard. | V3 |
| **V2** | Fix the hardcoded `region == "caudate"` in the 3 BA_only `02.pca.R` copies; add a guard asserting `nrow(pc.csv)` equals the region's design N. | V3 |
| **V6** | Add the `min(10^6, nrow())` guard to the DLPFC and hippocampus `02.pca.R`. | V3 |
| **V5** | Make `remove_ct_snps()` use its argument. | — |
| **A1** | Re-run `02` → `02b` → `02c` → `03` → `04` for all three regions, both branches. Archive the current `vmr.bed`, `top1_cpg.tsv` and `_m/pca/` first. | V1, V2, V3 |
| **A2** | Diff old vs new VMR sets per region: seed turnover, VMR-level turnover, boundary shifts, N per region. This is the number that determines how much downstream work is real. | A1 |
| **V4** | Decide X/Y policy and diagnose the caudate 3× sex-chromosome excess. | A2 |
| **V7** | Sort chunks numerically and assert FID/IID consistency in `02b`. | V3 |

### Stage B — re-verify everything downstream of VMR definition

| # | Task | Depends on |
|---|---|---|
| **B1** | Re-run the elastic-net predictability analysis on the corrected VMR sets (both branches). | A1 |
| **B2** | Re-run `meqtl-validation` modules 01–02 (CpG meQTL mapping, VMR burden) on corrected VMRs. | A1 |
| **B3** | Re-run module 10 (caudate downsampling) and module 13 (universe decomposition). Module 13's coverage-filter result is independent of all of this, but its sd-cutoff arm is not. | A1 |
| **B4** | Re-check the caudate gate. V2 is caudate-sparing and V1 inflates the bar, so the failed "not solely sample size" gate must be re-evaluated, not assumed to survive. | B2, B3 |
| **B5** | Re-derive the cross-region variability comparison, which module 13 could not reproduce (r = 0.377). Expect V1 + V2 to be the cause; confirm rather than assume. | A1, B3 |

### Stage C — previously open manuscript tasks (now mostly gated by Stage A/B)

| # | Task | Depends on |
|---|---|---|
| **#2** | Stop modules 04, 05, 08 propagating the failed caudate gate. The `02_summarize_phase4.py` fix is needed regardless of Stage A. | — (fix), B4 (verdict) |
| **#3** | Audit module 04 cross-region concordance denominators — direction concordance is exactly 1.000 ×3 with ORs 91/1237/22, which is degenerate. | B2 |
| **#4** | Re-derive the Phase 7 SCZ retain/supplement decision. | B2 |
| **#6** | Repackage supplementary data. | #2, #3, #4 |
| **#7** | Write the meQTL validation Results section. | #3, B4 |
| **#8** | Apply Experiment 1 positioning to manuscript `content/`. | — |
| **#9** | Rebuild main figures. | #3, #4, A2 |
| **#10** | Reconcile stale claim language. | #6–#9 |
| **#11** | Summarize meQTL results with `/analysis-to-pi`. | #2, #3, #4 |
| **#12** | Move status/planning markdown and the 19 per-module READMEs to the manuscript repo (needs `git rm` — merge `90d2456ba` tracked them). | #2–#5 |
| **#13** | Universe decomposition. **Done** — verdict `coverage_uniformity_not_sample_size`. Revisit only its sd arm, via B5. | done |
| **λ** | λ_GC = 3.884 in module 10 is still unaddressed. | — |

### Recommended order

V3 → V1 + V2 + V6 → A1 → **A2**. A2 is the decision point: it measures how much of the
existing downstream work actually changes, and everything in Stage B and C should be
scheduled against that number rather than assumed.
