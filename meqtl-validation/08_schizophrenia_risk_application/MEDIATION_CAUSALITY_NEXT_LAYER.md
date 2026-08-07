# Next layer: mediation and disease causality (Phase 7)

## Current evidence layer (v1)

```text
Level 1  Risk variant → CpG methylation          DONE (3 regions)
Level 2  VMR methylation ↔ expression/PSI        DONE (enrichment; not locus mediation)
Level 3  Same risk variant → expression/splicing DONE via GTEx (3/5 loci); LIBD internal FDR null
Level 4  Formal meQTL–eQTL/sQTL colocalization   NOT DONE (deferred)
Disease  Case–control molecular association      DONE (Analysis 7; secondary-VMR FDR only)
Causal   Mediation / MR                          NOT SUPPORTED (and not required)
```

v1 supports an **architecture / regulatory-application** claim. It does **not** support “methylation mediates schizophrenia risk” or “these variants cause disease via DNAm.”

---

## 1. Biological claim (what the next layer can honestly say)

**Realistic next claim (Level 3 + diagnosis):**

> At prioritized schizophrenia-risk loci, the same risk alleles that associate with CpG methylation in genetically anchored caudate VMRs also associate with nearby expression or splicing, and a subset of these molecular features differs by schizophrenia diagnosis after technical and exposure adjustment.

**Not yet defensible:** methylation mediates genotype→transcript or genotype→disease.

---

## 2. Primary hypothesis (falsifiable)

Among caudate SCZ-meQTL-supported VMRs with transcriptional coupling:

1. The index / credible-set risk variant (or strong LD proxy, r²≥0.8) is associated with the linked expression or PSI feature in the same donors or in a matched brain QTL resource.
2. Direction is biologically coherent with the methylation–transcript association.
3. (Secondary) VMR methylation and/or the linked transcript feature differs by schizophrenia diagnosis after covariates.

Fail if: no shared genetic support at any prioritized locus, or all diagnosis associations vanish after smoking/tox/antipsychotic adjustment.

---

## 3. Statistical unit and design

| Analysis | Unit | Design |
|---|---|---|
| Shared genetic support | SNP–feature pair within locus | Same AA donors if LIBD eQTL available; else GTEx/BrainSeq lookup |
| Colocalization | Locus (meQTL window ∩ eQTL/sQTL window) | Only loci with adequate regional signal |
| Diagnosis | Donor (VMR or feature) | Case–control among Phase 7 prioritized VMRs only |
| Mediation | Locus | Only after Level 3 + strong instrument; otherwise skip |

Primary region remains **caudate**. DLPFC/hippocampus are for concordance / specificity, not required for every locus.

---

## 4. Analyses that move one layer forward (ordered)

### A. Prioritize ≤5 illustrative loci (required gate)

From caudate FDR hits, rank loci by:

- fine-mapped / high-PIP risk variant (PGC3 FINEMAP / PIP tables)
- high VMR predictability
- internal meQTL FDR
- TX coupling (expression and/or PSI)
- mappability / non-segdup
- optional: external brain meQTL support

**Output:** `prioritized_loci.tsv` (≤5 rows) before any mediation language.

### B. Level 3 — shared genetic support (highest ROI)

For each prioritized locus, test whether the **same risk variant or r²≥0.8 proxy** associates with the linked gene/PSI:

1. **Internal (preferred if powered):** genotype → expression/PSI in LIBD AA caudate with the same covariate philosophy as meQTL.
2. **External (always):** lookup in GTEx v11 brain tissues  
   (`/projects/b1213/resources/public_data/gtex_v11/GTEx_Analysis_v11_eQTL/` and `.../sQTL/`, caudate / frontal cortex / hippocampus) and/or BrainSeq eQTL archive.
3. Record effect directions for risk allele on: methylation, expression/PSI, GWAS OR.

**Success:** ≥1 locus with risk→meQTL and risk→eQTL/sQTL (or strong LD proxy), coherent direction.  
This is the minimal step from “coupled” to “shared genetic support.”

### C. Analysis 5 — regional specificity (needed for main-text SCZ claim)

Already planned; now unblocked by multi-region Phase 7:

- Compare risk–CpG effect sizes / CIs across caudate, DLPFC, hippocampus
- Downsample caudate to matched N
- Shared-donor genotype×region where possible
- Do **not** call caudate-selective from significance pattern alone

### D. Analysis 7 — diagnosis validation (orthogonal disease layer)

For prioritized VMRs / linked features only:

```text
feature ~ primarydx + age + sex + ancestry PCs + technical + cell composition
```

Sensitivities: smoking, nicotine, cocaine, amphetamine, opioids, ethanol, antipsychotics, pH, PMI.

**Interpretation rule:** null case–control does **not** invalidate genetic meQTL (inherited effects exist in cases and controls). Positive diagnosis association is supportive, not causal proof.

### E. Formal colocalization (Level 4; conditional)

Run coloc/SuSiE-coloc **only** when:

- meQTL regional signal is strong
- eQTL/sQTL regional signal is strong
- LD reference matches ancestry reasonably
- single-signal assumption is plausible

Assets: PGC3 PIPs (`.../snps_in_gwas_loci.pip.gz`), GTEx SuSiE summaries, full PGC3 sumstats.  
In-repo coloc code is unfinished (`sex_context_brain` scaffolds) — expect new scripts, not reuse of a working pipeline.

**Do not** interpret prior-driven posteriors as disease mechanism.

### F. Mediation / MR (only if A–B succeed strongly)

If Level 3 is clear at ≥1 locus:

- Prefer **locus-level triangulation** over formal mediation in N≈150
- If mediation is attempted: two-step / product-of-coefficients with bootstrap; report as exploratory
- Weak-instrument MR and genome-wide mediation remain **non-goals**

---

## 5. Expected effect directions

| Link | Expected if architecture claim holds |
|---|---|
| Risk allele → methylation | Non-zero (already observed at FDR hits) |
| Higher predictability → higher SCZ-meQTL support | Positive (already observed) |
| Risk allele → linked expression/PSI | Non-zero; direction consistent with meth↔tx |
| Diagnosis → VMR/feature | Optional; may be null |

---

## 6. Minimal primary analysis for the next paper revision

1. Finish Phase 7 DLPFC + hippocampus + cross-region concordance table.  
2. Build ≤5 prioritized multi-omic loci.  
3. Level 3 shared genetic support (internal and/or GTEx).  
4. Analysis 7 for those loci only.  
5. Coloc only if powered; otherwise state Level 3 as ceiling.

**Main figure-worthy panel (if successful):** one locus track — GWAS / meQTL / VMR predictability / expression or PSI / cross-region effect forest.

---

## 7. Reviewer objections and planned responses

| Objection | Response |
|---|---|
| “This is mediation.” | No; Level 1–3 triangulation only unless formal mediation is powered. |
| “Caudate-only / N-driven.” | Analysis 5 downsample + cross-region effects. |
| “MethPC removed genetic signal / residual confounding.” | M3a locked with calibration; sensitivity without methPCs at prioritized loci. |
| “Wrong ancestry LD for coloc.” | Prefer internal LIBD QTL or ancestry-aware LD; down-weight European-only coloc. |
| “Case–control is drug/smoking.” | Prespecified tox/antipsychotic sensitivities. |
| “Index SNP ≠ causal.” | Use FINEMAP/PIP credible sets; treat index as tag. |

---

## 8. Go / no-go for stronger causal language

| Outcome | Language allowed |
|---|---|
| Level 3 at ≥1 prioritized locus + TX coupling | “Shared genetic support for methylation and transcription at SCZ-risk loci” |
| + Analysis 5 supports caudate biology | Main-text SCZ application |
| + Analysis 7 positive after sensitivities | “Also altered in cases” (still not mediation) |
| + Powered coloc | “Colocalized meQTL and eQTL/sQTL signals” |
| Formal mediation significant with strong instruments | Exploratory mediation sentence only |
| Fail Level 3 / Analysis 5 | Keep architecture claim; SCZ stays compact supplement |

---

## 9. Implementation order (do next)

1. Multi-region Phase 7 (submitted) → concordance summary script  
2. `06_prioritize_loci.py` (≤5 loci; FINEMAP/PIP + TX + predictability)  
3. `07_level3_shared_genetic_support.py` (GTEx lookup first; internal eQTL if available)  
4. `08_diagnosis_validation.py` (Analysis 7, prioritized only)  
5. Coloc module only after Level 3 shortlist is non-empty
