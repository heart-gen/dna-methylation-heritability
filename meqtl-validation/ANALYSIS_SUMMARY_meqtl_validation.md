## Analysis Summary: CpG cis-meQTL validation of VMR local genetic predictability

**Analysis unit:** `meqtl-validation/` (Phases 0–7 + Exp2/3 architecture depth; official TensorQTL downsample remap in `10_downsampling_caudate/` **complete**)  
**Role:** Consolidated technical and reproducibility summary for the VMR genetic-validation analyses  
**Status snapshot:** Phases 0–7 and the cross-region, donor-group, downsampling, and cell-type sensitivity analyses are complete  

### Purpose

Test whether variably methylated regions (VMRs) with greater local genetic predictability contain a greater burden of conventional CpG-level cis-meQTL evidence, and whether that genetically anchored architecture is reproducible across caudate, dorsolateral prefrontal cortex (DLPFC), and hippocampus and across Black American (AA) and white American (EA) donor groups—without treating BrainSEQ WGBS meQTL calls as independent external replication.

### Inputs

- Predefined VMRs and elastic-net local genetic-predictability summaries (AA discovery; EA portability from all-individuals elastic-net outputs).
- WGBS CpG methylation matrices and imputed genotypes for AA (primary) and EA (stratified) BrainSEQ donors.
- Locked Phase 1 covariates: **M3a** = age + sex + schizophrenia diagnosis + ancestry SNP PC1–5 + methylation PC1–5 (`config/covariates.yml`).
- Public brain meQTL resources catalogued in `inputs/data_dictionary/_m/public_meqtl_resources.tsv` (primary independent: Jaffe DLPFC 450K; Schulz hippocampus array).
- Existing VMR–expression and VMR–PSI association tables under `heritability/elastic_net_model/.../regulatory_context/`.
- Schizophrenia risk–locus inputs for Phase 7 application (`08_schizophrenia_risk_application/`).

### Methods Text

We mapped cis-meQTLs for CpGs within retained VMRs using TensorQTL cis permutation mapping with a ±500 kb window, MAF ≥ 0.05, and Storey *q*-value FDR controlled separately by brain region (`config/meqtl_parameters.yml`). The primary discovery cohort was Black American donors (caudate N = 153; DLPFC N = 111; hippocampus N = 116). After covariate-model sensitivity (M0–M5), we locked M3a as the primary AA model. For each VMR we aggregated CpG-level results into meQTL-burden metrics (notably the proportion of tested CpGs with FDR-significant cis-meQTL evidence) and modeled continuous local genetic predictability as the primary predictor, with technical/genomic adjustment and matched high- versus low-predictability analyses. Independent external support was tested using non-BrainSEQ public brain meQTL resources. Cross-region concordance, donor-group stratified mapping, transcriptional/splicing enrichment, and repeat/mappability sensitivity analyses were run as secondary architecture tests. For caudate sample-size sensitivity, we (i) retested full-N lead SNP–CpG pairs on 30 shared-donor-aware downsamples to N = 111 with BH-FDR (module `04`) and (ii) remapped the same sample lists with official TensorQTL permutation FDR (module `10`; **30/30 complete**). Shared-donor genotype × region models used within-region M3a residualization and cluster-robust SEs by donor.

### Results Text

**Internal CpG cis-meQTL discovery (AA, M3a).** TensorQTL tested 197,613 / 155,546 / 157,925 CpGs in caudate / DLPFC / hippocampus and identified 81,007 / 51,273 / 55,957 FDR-significant CpGs (FDR < 0.05). Genomic inflation λ_GC was high (≈4.98 / 3.72 / 3.92), consistent with abundant true cis signal in this design rather than a claim of well-calibrated null inflation alone.

**VMR meQTL burden validates predictability (primary claim).** In adjusted technical models, continuous local predictability associated positively with CpG meQTL burden in all three regions (coef ≈ 3.20 / 1.77 / 2.07; all *P* reported as 0 in model output tables, i.e. below floating-point resolution). Matched high- versus low-predictability contrasts remained positive (mean Δ proportion ≈ 0.82 / 0.61 / 0.69; permutation *P* ≈ 5.0 × 10⁻⁴ in each region).

**External validation.** Tissue-matched adjusted models supported the predictability gradient for Jaffe DLPFC (coef ≈ 0.149, *P* = 3.5 × 10⁻⁶) and Schulz hippocampus (coef ≈ 0.298, *P* = 2.4 × 10⁻¹⁷). BrainSEQ WGBS meQTL was excluded from independent external credit because of cohort overlap.

**Cross-region architecture.** Among shared tested CpGs, both-significant direction concordance was 0.912 / 0.905 / 0.924 and Pearson *z*-score correlations (either-significant) were 0.894 / 0.885 / 0.923 for caudate–DLPFC, caudate–hippocampus, and DLPFC–hippocampus. VMR meQTL-support Jaccard was modest (≈0.28–0.30 pairwise; 260 VMRs supported in all three). Higher predictability predicted shared VMR support (logistic OR ≈ 1.81 / 2.02 / 1.95; all *P* ≪ 0.05). Ninety-two donors had AA meQTL inclusion in all three regions.

**Caudate ≠ N (method-matched).** Official TensorQTL permutation-FDR remapping of caudate to N = 111 (30 reps; module `10_downsampling_caudate/`) yielded median **67,562.5** FDR-significant CpGs (IQR 67,181–67,922; median λ_GC ≈ 3.89), retaining a median **78.3%** of the 81,007 full-N FDR hits. All 30 replicates exceeded DLPFC (51,273) and hippocampus (55,957) Phase 1 n_sig (median ratios ≈ 1.32 / 1.21). Claim: caudate discovery excess is **not** explained solely by sample size (`tensorqtl_downsample_claim_snapshot.tsv`). Lead-SNP retention BH-FDR (module `04`; median retention 78.9%) remains a secondary sensitivity and must not be mixed with TensorQTL absolute n_sig. Shared-donor G×region screening of 3,000 stratified pairs found joint interaction FDR significance for 737 pairs (24.6%), indicating minority—not majority—regional effect heterogeneity.

**Donor-group (EA) depth.** EA stratified mapping (M0 covariates; N = 129 / 55 / 60) yielded 62,595 / 21,198 / 29,372 FDR CpGs. Predictability→burden remained positive (tech-adj coef ≈ 3.21 / 1.13 / 1.29). AA–EA predictability Spearman ρ ≈ 0.51 / 0.12 / 0.30. Among shared CpGs, both-significant direction concordance was 0.749 / 0.737 / 0.724. AA–EA discovery-rate gaps shrank after MAF/cis-SNP-density matching in DLPFC and hippocampus but not caudate, and were not eliminated; ancestry-specific biology is not claimed.

**Transcription / splicing.** meQTL-supported VMRs were enriched for expression associations in 3/3 regions (Fisher OR ≈ 1.88 / 12.8 / 8.18) and for PSI associations in 2/3 (OR ≈ 4.80 / 1.89; hippocampus not enriched, OR ≈ 0.61, *P* = 0.89).

**Repeat / repressive chromatin sensitivity.** Consolidated Phase 6 claims passed for H3K9me3 (3/3), quiescent chromatin (3/3), and predictability–burden (3/3); LINE/L1 enrichment was direction-consistent in 2/3 (DLPFC fragile under high-mappability filters).

**Schizophrenia-risk application (secondary).** In caudate AA M3a, 31 risk loci / 361 FDR risk–CpG pairs / 38 VMRs met targeted FDR; architecture enrichment toward higher predictability was positive (tech-adj coef ≈ 0.53, *P* = 1.6 × 10⁻⁴; matched Δ ≈ 0.096, perm. *P* ≈ 0.003).

### Figure and Table Notes

- Potential main figure: internal validation of predictability by CpG meQTL burden (Phase 2 adjusted + matched panels; regions as facets).
  - Rationale: primary genetic-validation claim.
  - Key message: continuous predictability predicts meQTL burden after technical adjustment and matching.
  - Required legend details: M3a covariates; FDR family = per region; matching variables.

- Potential main figure: cross-region concordance + donor-group portability (Phase 4/5).
  - Rationale: multi-region / diversity advance.
  - Key message: effect directions highly concordant; gradient portable; residual discovery gaps are not ancestry claims.
  - Include Exp3 method-matched TensorQTL downsample (module `10`: median n_sig 67,562.5; retention ≈0.78; exceeds DLPFC/hip).

- Potential main figure: repeat/repressive robustness (Phase 6 consolidated table).
  - Rationale: defends distal LINE/L1 / H3K9me3 architecture claim.
  - Note DLPFC LINE/L1 fragility.

- Potential supplementary figure/table: external Jaffe/Schulz models; EA M0 vs EA-M3a caudate sensitivity; lead-SNP retention (module `04`) vs official TensorQTL remap (module `10`); Phase 7 SCZ application if retained as proof-of-application.
  - Key columns for tables: region, n, coef, SE, *P*/FDR, model name, resource ID.

- Do not elevate BrainSEQ WGBS meQTL as independent external validation in figures.

### Reproducibility Information

- Analysis directory: `/projects/b1213/users/kynon/projects/dna-methylation-heritability/meqtl-validation/`
- Primary scripts: `01_cpg_meqtl_mapping/_h/04_tensorqtl_map.py`, `05_qc_summarize.py`; `02_vmr_meqtl_burden/_h/*`; `03_external_meqtl_validation/_h/*`; `04_cross_region_sharing/_h/*`; `05_donor_group_comparison/_h/*`; `06_transcription_splicing_integration/_h/*`; `07_repeat_mappability_sensitivity/_h/*`; `08_schizophrenia_risk_application/_h/*`; `10_downsampling_caudate/_h/*`
- Config: `config/meqtl_parameters.yml`, `config/covariates.yml`, `config/paths.yml`, `config/analysis_thresholds.yml`
- Key outputs:
  - AA M3a QC: `01_cpg_meqtl_mapping/{region}/_m/tensorqtl/qc/meqtl_qc_summary.tsv` (cis maps under `.../covariate_sensitivity/tensorqtl/M3a/`)
  - EA M0 QC: `01_cpg_meqtl_mapping/{region}/_m/tensorqtl/EA/qc/meqtl_qc_summary.tsv`
  - Burden: `02_vmr_meqtl_burden/_m/{region}/burden_model_results.tsv`, `matched_analysis_results.tsv`
  - External: `03_external_meqtl_validation/_m/phase3_criterion5_verdict.tsv`
  - Cross-region / Exp3: `04_cross_region_sharing/_m/phase4_claim_summary.tsv`, `caudate_downsample/downsample_claim_snapshot.tsv`, `gxregion/gxregion_claim_snapshot.tsv`
  - Official TensorQTL downsample: `10_downsampling_caudate/_m/tensorqtl_downsample_claim_snapshot.tsv`, `tensorqtl_downsample_replicate_results.tsv`, `tensorqtl_downsample_vs_regions.tsv`
  - Donor-group / Exp2: `05_donor_group_comparison/_m/aa_ea_effect_concordance_summary.tsv`, `maf_ld_matched_discovery_claim.tsv`, `experiment2_depth_claim_summary.tsv`
  - TX: `06_transcription_splicing_integration/_m/tx_enrichment_primary.tsv`, `phase5_claim_summary.tsv`
  - Repeat: `07_repeat_mappability_sensitivity/_m/phase6_claim_summary.tsv`, `consolidated_robustness_table.tsv`
  - SCZ: `08_schizophrenia_risk_application/_m/caudate/architecture_decision_snapshot.tsv`
- Log files inspected: Phase 1 TensorQTL SLURM logs under region `_m/logs/`; Exp3 downsample/G×region completion under `04_.../_m/logs/`; module `10` logs under `10_downsampling_caudate/_m/logs/`
- Execution command (representative): `sbatch` of module `step_*.sh` from each module `_h/` using conda env `/projects/p32505/opt/envs/genomics`
- Execution date: primary M3a lock and Phase 2–6 refresh 2026-08-01; Exp2/3 depth 2026-08-05; official TensorQTL remap **completed** (claim snapshot present; docs updated 2026-08-07)
- Git commit (repo HEAD at summary time): `b1c99ed36` (later uncommitted analysis outputs may exist beyond this commit)
- Workflow manager: SLURM batch scripts (no Snakemake/Nextflow for these modules)
- Compute environment: Northwestern Quest HPC; TensorQTL on `gengpu` A100; CPU steps on `genomics`
- Container or environment file: conda env `genomics` at `/projects/p32505/opt/envs/genomics` (no project-local `environment.yml` recorded for this module set)
- Python version: 3.11.13 (from `genomics` env query 2026-08-05)
- Key package versions: `tensorqtl` 1.0.10; `pandas` 2.3.3; `numpy` 2.2.6 (env query). Runtime logs did not consistently print package versions.
- Random seed: TensorQTL default seed 20260722 in Phase 1 map script; Exp3/10 downsample seed **20260805**
- Missing reproducibility information: per-job package-version dumps not standard in all SLURM logs; TensorQTL import warns that R/`qvalue` via `rfunc` is unavailable—code path uses `py_qvalue` in `04_tensorqtl_map.py` (confirm wording in Methods)

### Limitations and Integration Notes

1. High λ_GC in cis-meQTL maps is expected with dense true signal but should be described carefully; nonsignificant-site λ_NS informed M3a lock rather than raw λ_GC alone (`PHASE1_LOCK_DECISION.md`).
2. Elastic-net VMR scores are **local genetic predictability**, not calibrated locus heritability.
3. Lead-SNP retention BH-FDR and TensorQTL permutation FDR are different families—do not compare absolute n_sig across those methods in main text.
4. BrainSEQ WGBS meQTL must not be counted as independent Phase 3 support [@doi:10.1038/s41467-021-25517-3].
5. Donor-group differences in discovery rates are not labeled ancestry-specific; no formal genotype × donor-group interaction is claimed.
6. G×region interactions are significant for a minority of screened pairs; caudate discovery excess should not be framed as widespread region-selective effect heterogeneity from discovery counts alone.
7. LINE/L1 enrichment is fragile in DLPFC under high-mappability filters; cell-type composition of repeat-rich compartments remains unresolved in bulk WGBS.
8. Integrate the validated results into the manuscript with claims constrained by the documented sensitivity analyses and limitations.

### Claim checklist (analysis-side, current)

| Claim | Supported by outputs? |
|---|---|
| Predictability → CpG meQTL burden (AA, 3 regions) | Yes |
| Survives adjustment + matching | Yes |
| ≥1 independent public brain meQTL resource | Yes (Jaffe + Schulz) |
| Cross-region effect concordance | Yes |
| Donor-group gradient portability | Yes |
| Transcriptional coupling of meQTL-supported VMRs | Yes (expression 3/3; PSI 2/3) |
| Repeat/repressive enrichment with sensitivity | Yes (LINE/L1 2/3) |
| Caudate ≠ solely N (method-matched TensorQTL) | Yes (module `10`; median n_sig 67,562.5; retention 0.783) |
| Ancestry-specific biology | No (not claimed) |
