# Plan: Cell-type contribution to LINE/L1 / repressive-compartment enrichment

**Status:** Implemented — see `meqtl-validation/11_celltype_compartment_sensitivity/`  
**Decision:** `keep_main_figure_with_cell_composition_row` (`_m/celltype_claim_snapshot.tsv`)  
**Parent claim:** High local-genetic-predictability VMRs enrich in LINE/L1-rich and H3K9me3 / quiescent compartments (`07_repeat_mappability_sensitivity/`).  
**Why needed:** Bulk WGBS cannot distinguish true compartment biology from cell-composition shifts (especially oligodendrocyte / neuron fraction) that correlate with repeat-rich methylation.

---

## 1. Biological claim

We will test whether the association between continuous VMR local genetic predictability and LINE/L1 (and H3K9me3 / quiescent) overlap in caudate (primary; DLPFC/hippocampus secondary) remains after accounting for estimated cell-type composition, supporting the claim that genetically anchored methylation concentrates in repeat-rich repressive sequence beyond bulk cellularity.

## 2. Primary hypothesis

Among retained VMRs, LINE/L1 overlap (and H3K9me3 / quiescent labels) increases with continuous local genetic predictability after adjustment or matching that includes sample-level cell-composition covariates; the predictability → meQTL-burden relationship is unchanged.

**Fail if:** LINE/L1 enrichment vs predictability attenuates to null or reverses after cell-composition adjustment/matching in ≥2 regions, or is explained entirely by oligodendrocyte/neuron fraction proxies.

## 3. Statistical unit and design

| Layer | Unit | Design |
|---|---|---|
| Architecture enrichment | VMR | VMR-level annotation ~ z(predictability) + technical covariates + **sample-aggregated cell composition** |
| Sample-level confounding check | Donor/sample | Does cell fraction associate with mean methylation in LINE/L1 VMRs? |
| Orthogonal | CpG/VMR in external single-nucleus or sorted methylomes | Enrichment of high-predictability VMRs in cell-type-specific methylated compartments |

Primary region: **caudate** (largest N; strongest LINE/L1 signal historically).  
Secondary: DLPFC, hippocampus (note DLPFC LINE/L1 already fragile under high-mappability).

## 4. Primary comparison

**Predictor:** continuous VMR local genetic predictability (AA elastic-net score; same as Phase 6).  
**Outcomes (separate models):**
1. LINE/L1 overlap fraction (or binary LINE/L1-overlapping VMR)
2. H3K9me3 overlap / enrichment label
3. Quiescent chromatin label
4. (Control outcome) predictability → meQTL burden (must remain positive)

## 5. Covariates and cell-composition inputs

### Available now (in-repo)

| Source | Path | Use |
|---|---|---|
| RNA MuSiC proportions | `inputs/cell_proportions/_m/music-proportions-{region}.tsv` | Primary bulk RNA deconvolution (Astro, D1/D2-SPN, Immune, Inhib, Micro, OPC, Oligo) |
| MuSiC PCs (M5) | `01_cpg_meqtl_mapping/{region}/_m/covariate_sensitivity/covariates_M5.txt` | `cellPC1–3` already built for meQTL sensitivity |
| DNAm scMD proportions | `inputs/cell_proportions/_m/dnam-scmd-proportions-*.tsv` | Methylation-based deconvolution (n≈151 caudate AA) |
| DNAm cell PCs (M6d) | `.../covariates_M6d.txt` | `dnamCellPC1–3` |

### Exclusions / caveats

- Do **not** residualize VMR predictability on diagnosis or methPCs for this architecture test.
- Prefer composition **PCs** or a predeclared low-dimensional set (e.g., Oligo + neuron aggregate) to avoid overfitting with 8 collinear fractions.
- Document missingness (DNAm scMD drops ~2 caudate AA samples).

## 6. Expected effect direction

- Predictability positively associates with LINE/L1, H3K9me3, and quiescent labels after cell adjustment.
- Oligo fraction may correlate with global methylation; adjustment should shrink but not eliminate compartment enrichment if the claim is sequence/compartment-driven.
- meQTL-burden ~ predictability remains positive (sanity).

## 7. Minimal primary analysis

**Module:** `meqtl-validation/11_celltype_compartment_sensitivity/` (new)

1. Build VMR-level analysis table from Phase 6 `vmr_technical_annotations.tsv` + predictability + burden.
2. Attach **sample-mean** cell covariates? No — enrichment is VMR-level. Instead:
   - **Approach A (preferred):** re-estimate Phase 6 enrichment models with VMR methylation summaries that were residualized on cellPCs / dnamCellPCs at the sample level before VMR aggregation; compare coefficients to unadjusted.
   - **Approach B:** matched high vs low predictability VMRs stratified by genomic annotation, then test whether samples’ Oligo/neuron fractions differ between methylation quantiles of LINE/L1 VMRs (sample-level confounding).
   - **Approach C:** include region-level mean Oligo fraction as a cohort descriptor only (weak).

**Primary reportable test (Approach A):**  
For each region, compare consolidated robustness rows:

`LINE/L1 ~ predictability` before vs after cell-composition residualization of CpG/VMR methylation.

## 8. Sensitivity analyses

1. RNA MuSiC PCs vs DNAm scMD PCs (concordance of inference).
2. Adjust for Oligo fraction alone vs full PC set.
3. Restrict to high-mappability VMRs (already fragile in DLPFC).
4. Exclude segmental duplications / SNP-proximal CpGs (reuse Phase 6 filters).
5. Neuron-enriched vs glia-enriched annotation partitions if marker-gene / chromatin labels available.
6. Leave-one-cell-type-out: drop Oligo PC contribution and retest.

## 9. Orthogonal validation

Prioritize public **cell-type-resolved brain methylomes** (not new experiments):

1. **Libd / psychENCODE / Luo et al.–style snmC-seq or sorted NeuN± WGBS/array** with published cell-type DMRs or mean methylation tracks (hg38).
2. Test whether high-predictability VMRs overlap cell-type-specific hyper/hypomethylated regions more than low-predictability VMRs, **matched** on CpG density, length, and mappability.
3. If a caudate or striatum single-cell methylation atlas is available under `/projects/b1213/resources`, prefer that tissue; otherwise use cortical NeuN± as a limited external check and label as tissue-imperfect.

Do not claim cell-type of origin for a VMR without overlapping cell-resolved methylation or chromatin.

## 10. Main figure-worthy result

One consolidated panel/table row in the Phase 6 robustness table:

| Estimate | Original | +MuSiC cellPCs | +DNAm cellPCs | High-mappability | Direction consistent |
|---|---|---|---|---|---|
| LINE/L1 ~ predictability | … | … | … | … | Y/N |

Manuscript text only upgrades the compartment claim if ≥2 regions remain direction-consistent after cell adjustment.

## 11. Reviewer objections

| Objection | Response |
|---|---|
| “LINE/L1 = oligodendrocyte fraction” | CellPC / Oligo-adjusted models; sample-level Oligo vs LINE-VMR methylation check |
| “Deconvolution is noisy in AA bulk” | Dual RNA + DNAm deconvolution; report uncertainty; avoid single-fraction dogma |
| “Circular: methPCs include cell composition” | Architecture models use predictability from genotypes, not methPCs; cell adjustment is explicit |
| “DLPFC already fragile” | Prespecify caudate primary; DLPFC as sensitivity |
| “Need single-cell multiome” | Orthogonal public sn-methylation overlap; full multiome out of scope for the current analysis |

## 12. Go / no-go for manuscript language

| Outcome | Language |
|---|---|
| Enrichment survives MuSiC + DNAm cell adjustment in caudate + ≥1 other region | Keep LINE/L1–repressive claim in main Figure 5 with cell-composition row |
| Survives only in caudate | Main-text caudate; others supplemental |
| Attenuates to null after cell adjustment | Downgrade to “bulk compartment correlation; cell composition not excluded” |
| Opposite direction after adjustment | Remove causal/compartment wording; report as composition-sensitive |

## 13. Implementation order

1. **Audit** cell-proportion sample overlap with AA meQTL inclusion lists (write `11_.../_m/celltype_sample_overlap.tsv`).
2. **Residualize** prepared CpG bed (or VMR means) on M5 cellPCs and M6d dnamCellPCs; recompute VMR mean methylation.
3. **Rerun** Phase 6 enrichment + burden models on residualized matrices (`03_run_robustness_analyses.py` pattern).
4. **Append** rows to `consolidated_robustness_table.tsv` (`cellPC_adj`, `dnamCellPC_adj`).
5. **Sample-level check:** Oligo fraction ~ mean methylation of high-LINE VMRs.
6. **Optional orthogonal:** public cell-type methylation track overlap.
7. **Update** the consolidated robustness summary and technical documentation.

## 14. Explicit non-goals

- New single-cell library generation
- Cell-type-specific meQTL mapping genome-wide
- Claiming a VMR is “oligodendrocyte-specific” from bulk deconvolution alone
- Replacing mappability sensitivity with cell adjustment

## 15. Effort estimate

| Task | Effort |
|---|---|
| Overlap audit + residualized VMR means | 0.5–1 day |
| Phase 6 reruns + consolidated table | 1 day |
| Sample-level Oligo checks + writeup | 0.5 day |
| Public sn-methylation overlap (if assets exist) | 1–2 days |
| **Total to decision-ready** | **~3–4 days** |
