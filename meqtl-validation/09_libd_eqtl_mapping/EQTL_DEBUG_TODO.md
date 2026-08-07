# LIBD eQTL debug task list

**Status:** Open — genome-wide discovery is far below expectation (~1–2 eGenes at FDR 0.05; BrainSeq/GTEx caudate typically yield thousands).

**Module:** `meqtl-validation/09_libd_eqtl_mapping/`  
**Primary cohort lock (unchanged):** AA, Age > 13, Dx ∈ {Control, SCZD}, `TOPMed_LIBD.AA`, ±500 kb, MAF ≥ 0.01.

---

## Evidence already in hand

| Arm | Path | Genes | Expr PCs | eGenes FDR&lt;0.05 | λ_GC |
|---|---|---:|---:|---:|---:|
| CPM + filterByExpr | `_m/caudate/genes/` | 33,206 | 20 | 1 | 0.89 |
| RPKM mean&gt;0.2 | `_m/caudate/genes_rpkm/` | 18,473 | 16 | 2 | 0.89 |
| RPKM + max 5 PCs | `_m/caudate/genes_rpkm_pc5/` | 18,473 | 5 | 1 | 0.95 |

Comparison table: `_m/caudate/sensitivity_rpkm_vs_cpm.tsv`

**Ruled out (for FDR&lt;0.05 collapse):**
- [x] filterByExpr + TMM-CPM vs mean-RPKM filter/phenotype
- [x] Expression PC count 20 → 5 (no meaningful recovery)

**Concerning patterns:**
- Top leads are rare (AF ≈ 0.012–0.03; often ≤7 minor-allele carriers)
- Mild λ_GC deflation, not massive over-correction
- Level 3 still passes via **GTEx** external lookup; internal LIBD map is the weak link

---

## Debug tasks (ordered)

### A. Integrity / wiring (do first)

- [ ] **A1. Sample alignment audit**  
  Confirm phenotype BED columns, covariate `sample_id`, and `PgenReader.sample_ids` are the same BrNum set and order after `select_samples`.  
  Output: `_m/caudate/debug/sample_alignment.tsv` (n overlap, any permutation check).

- [ ] **A2. Dosage sanity on common SNPs**  
  Pick 5–10 common autosomal SNPs (MAF&gt;0.1); report mean dosage, variance, HWE-ish allele frequency vs `.pvar`/reference.  
  Fail if dosages are near-constant or AF wildly off.

- [ ] **A3. Chromosome / TSS matching**  
  Confirm phenotype `#chr` matches genotype chrom keys after `prepare_chr_matched_genotypes` (chr-prefix). Spot-check cis window SNP counts per gene on chr1/chr22.

- [ ] **A4. Phenotype matrix sanity**  
  Per-gene variance, fraction zero/constant, and that BED values match `normalized_expression.tsv.gz` for random genes/samples.

### B. Positive controls

- [ ] **B1. Curate ≥20 known caudate eGenes**  
  From GTEx v11 caudate eGenes and/or prior BrainSeq eQTL tables (prefer genes with common lead SNPs).  
  Output: `_m/caudate/debug/positive_control_egenes.tsv`

- [ ] **B2. Targeted SNP–gene tests**  
  For each control gene, test published/common lead SNP (or best local common SNP) with the same covariates as TensorQTL.  
  Success: majority show nominal association in expected direction.

- [ ] **B3. TensorQTL lookup for controls**  
  Pull `pval_nominal` / `pval_beta` / `qval` for control genes from existing cis outputs; quantify how many would be eGenes under looser thresholds.

### C. Filter / model sensitivities (after A–B)

- [ ] **C1. Zero expression-PC arm**  
  Covariates: `Sex + Dx + Age + snpPC1–3` only (no expr PCs). Remap cis on RPKM or CPM BED.  
  If still ~0 eGenes → PCs fully ruled out.

- [ ] **C2. MAF ≥ 0.05 remap**  
  Same phenotype/covariates as best arm from C1; raise TensorQTL `--maf 0.05`.  
  Goal: suppress rare-variant-dominated leads.

- [ ] **C3. Autosomes-only**  
  Drop chrX phenotypes (current RPKM FDR hits were X-linked); recompute eGene counts.

- [ ] **C4. Nominal QQ / histogram**  
  Genome-wide best-SNP `pval_nominal` and `pval_beta` QQ for one arm; save under `_m/caudate/debug/figures/`.

### D. External cross-checks

- [ ] **D1. Compare lead variants to GTEx caudate signif_pairs** for positive-control genes (variant identity / LD proxies).
- [ ] **D2. Optional:** rerun one chromosome with sex_context-style interaction-free standard model settings if a published LIBD FastQTL table exists for overlap.

### E. Decision / lock

- [ ] **E1. Root-cause note**  
  Write `_m/caudate/debug/ROOT_CAUSE.md` with go/no-go for using internal LIBD eQTL in Level 3 vs GTEx-only support.

- [ ] **E2. Update primary recipe**  
  If a sensitivity recovers normal eGene scale, promote that arm to `genes/` primary and refresh Level 3 (`08` step_6b).  
  If not, document internal map as failed QC and keep GTEx as Level 3 external support.

---

## Success criteria for “fixed”

At least one of:
1. ≥ hundreds of autosomal eGenes at FDR 0.05 with MAF≥0.05 leads and non-pathological λ_GC (≈0.9–1.1), **or**
2. Clear documented bugfix (alignment/dosage/chrom) with before/after eGene counts, **or**
3. Explicit decision that internal LIBD cis map is not usable and Level 3 relies on GTEx (+ targeted tests only).

---

## Out of scope for this debug

- Schizophrenia risk-variant architecture (stays in `08`)
- Expanding to EA/CAUC or multi-region until caudate AA map is credible
- Formal colocalization (Level 4) until internal or external eQTL signal is trustworthy at prioritized loci
