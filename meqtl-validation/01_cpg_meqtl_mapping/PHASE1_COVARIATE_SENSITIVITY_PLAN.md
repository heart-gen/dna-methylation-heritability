# Phase 1 covariate / latent-factor sensitivity — analysis plan

**Status:** plan only (not yet implemented)  
**Parent:** `PHASE1_LOCK_TODO.md`  
**Goal:** Choose and lock one primary CpG cis-meQTL covariate model before treating Phase 1–6 internal meQTL calls as final.

---

## 1. Biological claim

Residual inflation among non-significant CpGs (λ_NS ≈ 1.46–1.50 on `pval_beta`) may reflect incomplete technical confounding, mild polygenic/cis leakage, or true abundant signal. We will test whether a **prespecified, small** covariate/latent grid improves calibration **without** removing known genetic meQTL signal, then lock one primary model.

## 2. Primary hypothesis

**H1:** Adding a small number of methylation latent factors (and/or PMI/pH) reduces λ_NS toward ~1.0–1.1 while retaining (a) positive-control overlap with Jaffe/Schulz external meQTLs and (b) the Phase 2 predictability–burden relationship.

**H0 (lock current model):** No specification meaningfully improves λ_NS without material loss of genetic signal/external overlap; residual λ_NS is accepted as signal-driven and documented.

## 3. Statistical unit and design

| Element | Spec |
|---|---|
| Unit | Individual CpG (same phenotype BED as current Phase 1) |
| Samples | Same AA discovery inclusion lists (`preflight/sample_inclusion_primary.tsv`) |
| Regions | Caudate (pilot first), then DLPFC and hippocampus with the **same locked grid** |
| Cis window | ±500 kb (unchanged) |
| Engine | TensorQTL cis permutations + region-specific Storey qval |
| Phenotype | Non-residualized WGBS CpG methylation (do not residualize then re-test SNPs) |
| Genotypes | Existing filtered AA pfiles (unchanged) |

Do **not** change CpG/SNP inclusion mid-grid. Only the covariate matrix changes.

## 4. Primary comparison

Compare covariate models on the **same** CpGs and samples. Decision metric hierarchy:

1. λ_GC among `qval > 0.05` on `pval_beta` (**λ_NS**; primary calibration)
2. Retention of known external meQTL CpGs (Jaffe for DLPFC; Schulz for hippocampus; Jaffe overlap as secondary for caudate)
3. n significant at `qval ≤ 0.05` and lead-effect size distribution (signal retention)
4. Phase 2 burden ~ predictability direction/significance on the model’s significant set (caudate pilot)

## 5. Covariates and exclusions

### Current locked-candidate baseline (M0)

```text
agedeath + sex + primarydx + snpPC1–5
```

N: caudate 153 / DLPFC 111 / hippocampus 116. No PMI, pH, batch, cell composition, or latent factors today.

### Prespecified model grid (small; no large SV dump)

| ID | Specification | Rationale |
|---|---|---|
| **M0** | Baseline (current) | Reference |
| **M1** | M0 + `pmi` + `ph` | Observed technical covariates (near-complete in AA) |
| **M2** | M0 + snpPC6–10 | Extra ancestry structure |
| **M3a** | M0 + PEER/k=5 | Small latent set |
| **M3b** | M0 + PEER/k=10 | Mid latent set |
| **M3c** | M0 + PEER/k=15 | Upper bound; stop here |
| **M4** | M1 + PEER/k\* | Best technical + best k from M3\* (chosen after M3, before locking) |
| **M5** | M0 + top cell-proportion PCs (≤3) | Sensitivity only (`covariates.yml`: cell composition = sensitivity) |

**Method note:** The genomics env lacks the `peer` package; latent factors are
**PCA on M0-residualized CpG phenotypes** (`methPC1–k`), documented as the
planned PEER-style substitute in `latent_estimation_summary.tsv`.

**Exclusions from this lock decision:** toxicology/smoking grids, diagnosis interactions, EA expansion — those are separate sensitivity tracks, not Phase 1 lock criteria.

**Missing data rule:** If `pmi`/`ph` missingness would drop >5% of samples, impute median within region×sex for the sensitivity matrix and record; do not silently drop samples relative to M0 unless unavoidable.

## 6. Expected effect direction

- λ_NS decreases as latent/tech factors are added (if confounding drives inflation).
- External overlap (OR or enrichment of published meQTL CpGs among internal significant CpGs) stays stable or improves.
- n significant may drop modestly; a **large** drop with flat external enrichment suggests over-correction.
- Phase 2 predictability → burden should remain positive in ≥2 regions if the model is valid.

## 7. Minimal primary analysis

1. **Pilot region = caudate** (largest N; strongest current λ).
2. Build covariate matrices `covariates_{model_id}.txt` for M0–M3c (+ M1).
3. Rerun TensorQTL cis for each model (GPU preferred; same chunking/FDR).
4. For each model compute:
   - λ_GC overall and λ_NS (`qval > 0.05`)
   - QQ (`pval_beta`)
   - n sig (`qval ≤ 0.05`)
   - overlap enrichment vs Jaffe (and Schulz) CpGs present in the VMR CpG set
5. Select candidate k / model using caudate metrics.
6. Confirm **DLPFC** and **hippocampus** with M0 vs chosen candidate (+ M1 if not already nested).
7. Write lock decision memo; update `config/covariates.yml` primary formula.

## 8. Sensitivity analyses (within this lock package)

- M5 cell-composition (report; do not prefer over PEER unless clearly superior on λ_NS **and** external retention).
- M2 snpPC6–10 (ancestry vs methylation latent contrast).
- Optional: high-mappability CpG subset λ_NS for the locked model only (ties to Phase 6; not required to choose the model).

## 9. Orthogonal validation

- **External CpG overlap:** Jaffe (primary for DLPFC; secondary elsewhere), Schulz (primary for hippocampus).
- **Architecture retention:** re-run Phase 2 adjusted_minimal burden model on the locked significance calls for all three regions; require positive predictability coef in ≥2 regions.
- BrainSeq full catalogs remain **out of scope** (cohort overlap; not independent).

## 10. Main figure-worthy result

One comparison panel/table:

- x: model ID  
- y: λ_NS  
- annotations: n sig, external enrichment OR, Phase 2 burden coef  

Companion QQ facets for M0 vs locked model.

## 11. Reviewer objections and responses

| Objection | Response |
|---|---|
| λ_all ≫ 1 means the model is broken | Stratify: λ_NS ≈ 1.5; λ among significant is expected to be large with abundant true cis signal |
| PEER removes genetic signal | Require external enrichment retention and Phase 2 direction retention before accepting PEER |
| Cell composition is the real confounder | M5 is sensitivity-only; bulk meQTL literature often prefers PEER/SVs over noisy deconvolution |
| Why not many SVs? | Prespecified k∈{5,10,15} avoids automatic inclusion of a large SV set |
| Platform mismatch for external overlap | Enrichment uses CpG presence in published catalogs, not effect-size equality; tissue-matched primary tests prioritized |
| Changing calls invalidates Phase 2/3 | Recompute burden (and note Phase 3 external test is mostly VMR-predictability–based and should be stable); document deltas |

## 12. Implementation notes

### Suggested module layout

```text
meqtl-validation/01_cpg_meqtl_mapping/
  _h/
    09_estimate_latent_factors.py      # PEER/PCA → prepared/latent/
    10_prepare_covariate_models.py     # write covariates_{model}.txt
    11_compare_covariate_models.py     # λ_NS, n_sig, external overlap table
    12_lock_decision.py                # write PHASE1_LOCK_DECISION.md + tsv
    step_8_latent.sh
    step_9_covar_models.sh             # array over model×region (or caudate-first)
    step_10_compare_models.sh
  _m/{region}/covariate_sensitivity/
    covariates_M0.txt … covariates_M4.txt
    tensorqtl/{model}/…
    comparison_summary.tsv
  PHASE1_LOCK_DECISION.md              # created at end
```

### Compute strategy

- **Do not** overwrite current `tensorqtl/` primary outputs until a model is locked.
- Write sensitivity TensorQTL under `covariate_sensitivity/{model}/`.
- Caudate pilot first (3–5 GPU jobs for M0–M3c); then 2 regions × 2 models (M0 + winner).
- Reuse existing phenotype BEDs and genotypes; only `--covariates` path changes in `04_tensorqtl_map.py` (add CLI override if missing).

### Go / no-go for locking

**Lock new model (e.g. M1 or M3\*/M4) if:**

1. λ_NS improves by ≥0.15 **or** reaches ≤1.20, **and**
2. External enrichment OR does not drop >20% relative to M0 (tissue-matched resource), **and**
3. Phase 2 predictability coef remains >0 with p<0.05 in ≥2 regions.

**Lock M0 (current) if:**

- No model meets (1)–(3), **and** λ_NS remains ~1.45–1.55 with strong external/Phase 2 support — document residual inflation as likely true-signal burden rather than unmodeled batch.

**Do not lock / escalate if:**

- Latents collapse n sig by >50% with loss of external enrichment, **or**
- Phase 2 burden validation disappears under the candidate calls.

### Out of scope for this plan

- Re-deriving VMRs; changing FDR family; ±1 Mb cis window; EA-stratified meQTL; toxicology-adjusted meQTL; BrainSeq as external truth.

---

## Decision log (fill after runs)

| Date | Region | Chosen model | λ_NS | n sig | Ext. OR vs M0 | Phase 2 retained? | Decision |
|---|---|---|---:|---:|---:|---|---|
| | caudate | | | | | | |
| | dlpfc | | | | | | |
| | hippocampus | | | | | | |
