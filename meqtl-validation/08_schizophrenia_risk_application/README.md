# Phase 7: Schizophrenia-risk locus application

Biology-first manuscript summary: [`ANALYSIS_SUMMARY.md`](ANALYSIS_SUMMARY.md).

## Status

**Complete.** This is a focused application of the VMR genetic-predictability framework and is not required for the core Phases 1–6 validation.

Primary region: **caudate** (AA, locked M3a covariates).  
Analyses **1–7** + Level 3 (GTEx) are done. Level 4 coloc and formal mediation remain deferred.

**Retain decision:** `retain_main_text_proof_of_application` — see [`_m/decision/PHASE7_DECISION.md`](_m/decision/PHASE7_DECISION.md).

## Inputs (prespecified; no methylation-driven locus definition)

| Resource | Path |
|---|---|
| PGC3 locus intervals (hg38) | `/projects/b1213/resources/gwas/pgc3/fine_mapped_loci/to_hg38/_m/gwas_loci_hg38.tsv` |
| PGC3 index SNPs (Nature Supp Table 1) | `_m/raw/pgc3_primary_gwas_index_snps.tsv` (from `primaryGWAS_indexSNP.xls`) |
| AA M3a covariates / phenotypes / genotypes | `01_cpg_meqtl_mapping/caudate/_m/prepared/` + `genotypes/meqtl_AA` |
| Genome-wide cis leads (for VMR linking) | `01_.../caudate/_m/tensorqtl/qc/lead_snp_per_cpg.tsv.gz` |
| VMR predictability / burden | `02_vmr_meqtl_burden/_m/caudate/vmr_meqtl_burden.tsv.gz` |

## Run

```bash
cd meqtl-validation/08_schizophrenia_risk_application/_m
mkdir -p logs
# Per region (caudate / dlpfc / hippocampus):
J1=$(sbatch --parsable --export=ALL,REGION=caudate ../_h/step_1.sh)   # define loci + link VMRs
J2=$(sbatch --parsable --dependency=afterok:${J1} --export=ALL,REGION=caudate ../_h/step_2.sh)  # targeted meQTL + architecture + tx
# After all three regions finish:
sbatch ../_h/step_3.sh   # cross-region concordance + ≤5 locus prioritization
sbatch ../_h/step_4.sh   # Tier A: caudate donor-downsample targeted meQTL (N=111, 30 reps)
sbatch ../_h/step_5.sh   # Analysis 5: shared-donor genotype × region
# After prioritization + Level 3:
sbatch ../_h/step_7.sh   # Analysis 7 diagnosis + retain/supplement decision
sbatch ../_h/step_8_locus_panels.sh   # tidy locus tables + manuscript panels (rnaseq R)
```

## Scripts

| Step | Script | Analysis |
|---|---|---|
| 1 | `01_define_scz_loci.py` | Prespecify loci; match index SNPs to genotype panel (rsID → hg38) |
| 1 | `02_link_vmrs_to_loci.py` | Link VMRs via lead-meQTL / window / overlap; proximity = exploratory |
| 2 | `03_test_risk_variant_cpg_meqtl.py` | Targeted risk-variant–CpG tests; **separate FDR family** |
| 2 | `04_architecture_predictability.py` | Risk-meQTL support ~ continuous predictability (+ matched perm) |
| 2 | `05_tx_integration.py` | Enrichment of SCZ meQTL VMRs among expression/PSI links |
| 3 | `06_cross_region_concordance.py` | Effect concordance, locus sharing, √N discovery rates |
| 3 | `07_prioritize_loci.py` | Rank ≤5 illustrative loci (FINEMAP + meQTL + TX + mappability) |
| 4 | `08_make_caudate_downsample_lists.py` | Shared3-aware caudate sample lists (N=111) |
| 4 | `09_run_caudate_downsample_meqtl.py` | Retest same risk–CpG pairs on each downsample |
| 4 | `10_summarize_caudate_downsample.py` | Compare to DLPFC/hippocampus; claim snapshot |
| 5 | `11_shared_donor_gxregion.py` | Shared-donor G×region interactions (clustered SE) |
| 7 | `19_diagnosis_validation.py` | Analysis 7: Dx ~ prioritized VMR means + linked expr/PSI |
| 7 | `20_phase7_decision.py` | Consolidate §12.12–12.14 retain/supplement decision |
| 8 | `21_prepare_locus_panel_data.py` | Assemble GWAS / meQTL / forest / VMR / TX tidy tables |
| 8 | `22_plot_locus_panels.R` | Manuscript panels A–D (PDF+PNG) |

**Locus panels:** `_m/locus_panels/figures/` — heroes **rs8048039**, **rs13331198**; others supplemental. Panel key in `_m/locus_panels/figures/LOCUS_PANEL_README.md`.

## Multiple testing

- Genome-wide Phase 1 FDR family: unchanged
- Phase 7 risk-variant–CpG family: `significance.scz_risk_variant_cpg_fdr` (default 0.05) in `config/analysis_thresholds.yml`

## Prespecified decision rules

| Outcome | Action |
|---|---|
| Meets §12.12 primary success criteria | Retain in main manuscript |
| Conditional success (§12.13) | Proof-of-application / supplement |
| Stop criteria (§12.14) | Move to supplement or omit |

See `architecture_decision_snapshot.tsv` and `tx_decision_snapshot.tsv` after step_2.

## Explicit non-goals

- Schizophrenia PRS
- Another broad S-LDSC screen
- Nominal DMR × GWAS overlap as primary evidence
- Genome-wide genotype × diagnosis interactions
- Weak-instrument mediation / underpowered colocalization

## Cross-region v1 (AA M3a)

| Region | FDR pairs | Loci | VMRs | Pred. enrichment | TX-coupled VMRs |
|---|---:|---:|---:|---|---:|
| Caudate | 361 | 31 | 38 | Pass (matched too) | 8 |
| DLPFC | 217 | 17 | 18 | Pass (regression; matched null) | 0 |
| Hippocampus | 196 | 16 | 16 | Fail (tech-adj ns) | 0 |

See `_m/phase7_cross_region_summary.tsv`.

### Cross-region concordance + prioritization (`step_3`)

| Output | Path |
|---|---|
| Pairwise concordance | `_m/cross_region/pairwise_concordance.tsv` |
| Locus × region matrix | `_m/cross_region/locus_by_region.tsv.gz` |
| Claim snapshot | `_m/cross_region/cross_region_claim_snapshot.tsv` |
| Top 5 loci | `_m/prioritized/prioritized_loci.tsv` |
| Rationale | `_m/prioritized/prioritized_loci_rationale.tsv` |

### Tier A caudate downsample (`step_4`)

Shared-donor-aware downsampling of caudate AA to N=111 (30 replicates): keep 92 all-region donors, fill randomly, retest the full caudate risk–CpG pair family with M3a covariates.

| Output | Path |
|---|---|
| Design / lists | `_m/caudate_downsample/downsample_design_summary.tsv` |
| Per-rep discovery | `_m/caudate_downsample/downsample_replicate_results.tsv` |
| vs DLPFC/hip | `_m/caudate_downsample/downsample_vs_regions.tsv` |
| Claim snapshot | `_m/caudate_downsample/downsample_claim_snapshot.tsv` |
| Prioritized stability | `_m/caudate_downsample/downsample_prioritized_stability.tsv` |

### Shared-donor genotype × region (`step_5`)

For donors with all three regions, residualize within region on M3a, then fit  
`resid ~ g + C(region) + g:C(region)` with donor-clustered SEs (caudate reference).

| Output | Path |
|---|---|
| Pair results | `_m/gxregion/gxregion_pair_results.tsv.gz` |
| Summary | `_m/gxregion/gxregion_summary.tsv` |
| Prioritized rollup | `_m/gxregion/gxregion_prioritized_locus_summary.tsv` |
| Claim snapshot | `_m/gxregion/gxregion_claim_snapshot.tsv` |

Genome-wide caudate TensorQTL downsample (method-matched perm-FDR) is complete in [`../10_downsampling_caudate/`](../10_downsampling_caudate/) (30× N=111; Pass).

## Level 3 — shared genetic support (`step_6*`)

Risk variant → expression for prioritized loci. Genome-wide LIBD eQTL mapping lives in
[`../09_libd_eqtl_mapping/`](../09_libd_eqtl_mapping/); this module only builds targets,
looks up GTEx, and tests index SNPs against linked genes.

```bash
# 1) Genome-wide LIBD AA eQTL (reusable module)
cd meqtl-validation/09_libd_eqtl_mapping/_m
J9=$(sbatch --parsable --export=ALL,REGION=caudate ../_h/step_1.sh)
sbatch --dependency=afterok:${J9} --export=ALL,REGION=caudate ../_h/step_2.sh

# 2) SCZ Level-3 consumers
cd ../../08_schizophrenia_risk_application/_m
J6A=$(sbatch --parsable ../_h/step_6a_level3_targets_gtex.sh)
sbatch --dependency=afterok:${J6A} ../_h/step_6b_level3_eqtl_tests.sh
```

| Output | Path |
|---|---|
| Gene targets | `_m/level3/level3_gene_targets.tsv` |
| GTEx lookup | `_m/level3/gtex/` |
| LIBD cis map (canonical) | `../09_libd_eqtl_mapping/_m/caudate/genes/` |
| Risk-variant eQTL | `_m/level3/libd_risk_variant_eqtl.tsv.gz` |
| Claim snapshot | `_m/level3/level3_claim_snapshot.tsv` |

## Next layer (mediation / causality)

See [`MEDIATION_CAUSALITY_NEXT_LAYER.md`](MEDIATION_CAUSALITY_NEXT_LAYER.md) for the evidence ladder and required analyses to move beyond architecture claims.

### Analysis 7 — diagnosis (`step_7`)

Prioritized loci only (caudate AA). Primary model:
`feature ~ primarydx + agedeath + sex + snpPC1–5`
(+ cellPC / tox / PMI / pH sensitivities).

| Output | Path |
|---|---|
| Associations | `_m/diagnosis/diagnosis_association_results.tsv` |
| Locus rollup | `_m/diagnosis/diagnosis_locus_rollup.tsv` |
| Claim snapshot | `_m/diagnosis/diagnosis_claim_snapshot.tsv` |
| Retain decision | `_m/decision/PHASE7_DECISION.md` |

Result snapshot: best-VMR FDR null; secondary TX-linked VMR 9804 (locus rs8048039) FDR-significant and survives smoking / lifetime antipsychotic nominal sensitivity; linked transcript Dx associations null. Null/weak Dx does not invalidate genetic meQTL evidence.

## Deferred beyond current package

- Formal colocalization (Level 4)
- Formal mediation / MR
- LD-proxy R² computation beyond index-SNP identity (PGC LD-friends retained in raw table for follow-up)
- LIBD genome-wide eQTL QC repair (`09_libd_eqtl_mapping/EQTL_DEBUG_TODO.md`) — Level 3 currently via GTEx
