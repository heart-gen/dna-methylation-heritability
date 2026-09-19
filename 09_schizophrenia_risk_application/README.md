# 09_schizophrenia_risk_application — required disease application

Tests whether schizophrenia-risk variants regulate methylation within genetically anchored VMRs. Intended for the main text, conditional on surviving corrected VMRs and the new local-genetic-control axis.

**Status: implemented, not accepted.** Upstream gates are satisfied — `05_cpg_meqtl_burden` (`cmb-AA-*-20260825`) and `07_transcription_splicing_coupling` (`tsc-AA-*-20260902`) both record passing acceptance gates — so the module is unblocked (AGENTS.md §6). No run of this module has been accepted; nothing here may be cited.

## Migrating from

`meqtl-validation/08_schizophrenia_risk_application/`, whose Phase 7 decision file records `retain_main_text_proof_of_application`.

See `MIGRATION_MANIFEST.tsv` for the legacy paths, their downstream consumers,
and retirement status. Legacy directories stay in place until their row reads
`validated_replacement`.

## That decision does not carry forward

The existing Phase 7 result (31 caudate loci, 361 pairs, 38 VMRs, eight
TX-coupled VMRs) was conditioned on the legacy predictability metric and
pre-repair VMR sets. Per AGENTS.md §8 it is a hypothesis to retest. Hero loci
`rs8048039` and `rs13331198` remain **candidates**: the prioritization rule in
`config/schizophrenia.yml` names no locus, and `_h/10_prioritize_loci.R`
contains none.

## Preserved design principles

PGC schizophrenia loci are defined independently of methylation results —
`_h/01_define_scz_loci.R` reads only the published intervals, the published
index-SNP table and the public summary statistics, and records
`methylation_used = FALSE`. Risk-variant×CpG tests keep their own FDR family,
corrected once over all autosomes in `_h/03b_combine_risk_variant_tests.R`.
Association is tested against the relative local SNP contribution score, not
absolute PVE or legacy predictability. GTEx eQTL evidence is support, not proof
of mediation. At most five illustrative loci are prioritized by a prespecified
rule.

The integration analysis asks whether schizophrenia-linked VMRs are enriched
for LINE/L1, H3K9me3, quiescent chromatin, high-mappability repeat intervals,
and expression/splicing coupling. Positive connects the disease application to
the repeat/repressive architecture; null presents Phase 7 as a separate proof of
disease relevance. Both are results.

## Colocalization

`config/analysis_thresholds.yml:phase7_scz` listed coloc as `deferred` while no
ancestry-matched QTL resource was wired up. PI decision 2026-09-09 authorised
it, and AGENTS.md §7.8 permits the claim once the analysis "has been run with
adequate ancestry-matched LD and passes its own gate". Two arms:

| Arm | Pair | LD | Status |
| --- | --- | --- | --- |
| `gtex_eqtl`, `gtex_sqtl` | PGC3 European GWAS × GTEx v11 European brain | matched | gate-bearing; may support a claim |
| `meqtl` | PGC3 European GWAS × this cohort's African-American CpG meQTL | **not matched** | exploratory only; never claimable |

The meQTL arm is the biologically interesting pair, but coloc assumes both
studies share an LD structure and these do not. Every row of that arm carries
`ld_ancestry_matched = FALSE` and `arm_status = CROSS_ANCESTRY_LD_UNMATCHED`,
is excluded from the gate, and cannot set `claimable_colocalization` — a
condition the gate re-checks rather than trusting. It becomes gate-eligible
without a code change once the `all_individuals` cohort carries EA donors across
all three regions: set `colocalization.arms.meqtl.qtl_ancestry: EUR` and
`gate_eligible: true`.

`coloc.abf` is primary; `coloc.susie` runs as a sensitivity check on the
single-causal-variant assumption over prioritized loci only, and never
overturns the primary result.

## Pipeline

| Stage | Script | Purpose |
| --- | --- | --- |
| 00 | `00_new_run.R` | Open the run; gate upstreams 01/02/04/05/06/07; assert one shared `vmr_set_id`; freeze the coloc arms |
| 01 | `01_define_scz_loci.R` | PGC3 intervals + index SNPs, lifted to hg38; GWAS sliced to locus windows (via `slice_gwas_sumstats.sh`) |
| 02 | `02_link_loci_to_vmrs.R` | Link loci to corrected VMRs; freeze the tested universe and the background set |
| 03 | `03_risk_variant_cpg_meqtl.py` | Per autosome: index SNPs + LD proxies × member CpGs, from Module 05's nominal pairs |
| 03b | `03b_combine_risk_variant_tests.R` | Pool the autosomes; apply the module's own FDR family once |
| 04 | `04_architecture_axis.R` | Enrichment of SCZ-linked VMRs along `local_snp_contribution_score_z` |
| 05 | `05_integration_enrichment.R` | The five prespecified repeat/repressive annotations, with per-region constraints |
| 06 | `06_transcriptional_coupling.R` | Project Module 07's accepted coupling onto loci |
| 07 | `07_gtex_support.py` | External shared-variant support from GTEx v11 brain eQTL/sQTL |
| 08 | `08_prepare_coloc_regions.py` | Per autosome: harmonise GWAS and QTL statistics onto common variants |
| 09 | `09_run_coloc.R` | Per autosome: `coloc.abf` per region, per arm |
| 09b | `09b_combine_coloc.R` | Pool coloc; roll up to loci, gate-eligible arms kept separate |
| 10 | `10_prioritize_loci.R` | Apply `ranked_composite_v1`; ≤5 loci, every ranking key recorded |
| 11 | `11_coloc_susie_sensitivity.R` | `coloc.susie` credible-set sensitivity on prioritized loci |
| 12 | `12_apply_gates.R` | Conduct gate + coloc gate; retention criteria; interpretation constraints |
| 13 | `13_plot_locus_panels.py` | Stacked regional panels for the prioritized loci |
| 14 | `14_finalize_run.R` | Seal: session info, manifest fields, checksums, read-only |
| 15 | `15_cross_region_axis_concordance.R` | **Module-level, not per-run.** Cross-region axis concordance over the three accepted per-region runs; resolves decision 2; measures donor overlap and upstream freshness. Writes `_m/combined/` |
| 16 | `16_axis_downsampling_sensitivity.R` | **Module-level.** Repeats the axis contrast in the three `AA.n118r*` caudate cells: does the contrast survive matching n to 118? Writes `_m/combined/` |

## Two independent decisions

PI 2026-09-18. This module's framing question is *does regional variation in
local genetic control of methylation intersect schizophrenia-relevant regulatory
biology, and is that relationship shared or region-dependent?* That has two
separable answers, and Module 08 tier 3 speaks to only one of them. They are
recorded as separate fields in `scz-decision.tsv` and
`_m/combined/scz-application-decisions-{cohort}.tsv`:

| decision | question | gated on | scope |
|---|---|---|---|
| `caudate_magnitude_claim` | Is caudate's *stronger* signal more than its larger donor count? | Module 08 tier 3 reading, and nothing else | caudate runs only |
| `scz_application_retention` | Is there a disorder-related pattern across all three regions? | H1 concordance (stage 15) + the per-region criteria, **excluding** `caudate_not_sample_size_artifact` | module-level |

**H1, as locked:** SCZ-linked VMRs show *lower* local genetic-control scores in
all three regions, with statistical support in at least two, including at least
one non-caudate region. A negative-but-nonsignificant third region does **not**
falsify H1.

Three constraints that the code enforces and the manuscript must respect:

- **The three regions are not independent replicates.** 100 donors appear in all
  three; DLPFC and hippocampus share 115 of 118 (Jaccard 0.96); the union is 168
  donors. Stage 15 computes this and writes it out. Report *concordance*, never a
  pooled p-value — the fixed-effect estimate stage 15 emits is labelled
  descriptive-only and its SE is anticonservative by construction.
- **Caudate may contribute to concordance but never carry it.** It is perfectly
  confounded with sequencing batch (AGENTS.md §8.1), hence
  `axis_requires_support_outside: caudate`.
- **Attenuation is not bias.** Module 08 tier 3 shows donor count is a plausible
  major contributor to caudate's larger *magnitude*. It does **not** show the
  caudate estimate is biased. Write "caudate magnitude attenuates after
  n-matching", never "n-inflated".
- **No blanket repressive-chromatin claim.** The direct linked-vs-background
  annotation test exists (stage 05), but claimability is narrower than
  significance. Stage 15 reports, per annotation, the regions where the depletion
  is claimable; cite those, and do not generalize to "repressive chromatin".

Submit one cell:

```
./_h/submit_schizophrenia.sh AA caudate      # DRY_RUN=1 to print the job graph
```

Stages 01–02 run on the submit host; 03 and 08–09 are 22-way arrays.
Colocalization runs in the `coloc` conda env (`00_shared/slurm.sh:run_r_coloc`);
`coloc` and `arrow` are not in `epigenomics`.

## Acceptance gate

`_h/12_apply_gates.R` writes `results/scz-decision.tsv`. It checks that the
analysis was **conducted** correctly, not that it produced a positive result: a
null Phase 7 is a reportable finding (AGENTS.md §7.8 provides for it
explicitly), and Modules 06 and 07 set the precedent that a reportable null
seals.

1. Locus definition recorded `methylation_used = FALSE`.
2. Every stage produced its output table.
3. The risk-variant FDR family is the one config names, corrected once.
4. At least `min_loci_tested` loci, `min_vmrs_linked` VMRs and
   `min_pairs_tested` pairs entered the analysis.
5. At least one architecture test produced a finite estimate.
6. Colocalization ran on at least `min_loci_coloc_tested` ancestry-matched loci
   with at least `min_variants_shared` harmonised variants each.
7. No claimable colocalization originates from an ancestry-unmatched arm.

Decision codes: `PASS_SCZ_APPLICATION_QC`, `PASS_SMOKE_ONLY_NOT_ACCEPTABLE`,
`FAIL_SCZ_APPLICATION_QC:<reasons>`.

### Main-text retention is reported separately

The five prespecified criteria in `config/schizophrenia.yml:retention_criteria`
are evaluated into `results/retention-criteria.tsv`, each `PASS`, `FAIL`,
`NOT_APPLICABLE_NON_PRIMARY_REGION`, or `PENDING_MODULE_08`.

**`caudate_not_sample_size_artifact` is always `PENDING_MODULE_08`.** The
caudate downsampling arm belongs to `08_region_donor_generalization`, which is
not implemented, so the criterion is *unevaluable* — neither satisfied nor
failed. `main_text_retention` therefore stays `PENDING_MODULE_08` regardless of
how the other four resolve, so the open dependency cannot be lost in the
writing. When Module 08 lands, set `gates.require_module_08_downsampling: true`
and wire its downsampling result into `_h/12_apply_gates.R`.

Note also that Module 06's accepted S-LDSC result is **null**
(`sldsc_supports_brain_enrichment = FALSE`), which
`config/analysis_thresholds.yml` lists as an `omit_or_supplement_if` condition.
The decision file carries that upstream value.

## Accepted runs

| run_id | cohort | region | vmr_set_id | accepted_on | accepted_by | decision | notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| scz-AA-caudate-20260918 | AA | caudate | vmrset-AA-caudate-937a41979978 | 2026-09-19 | Kynon J. Benjamin | PASS_SCZ_APPLICATION_QC | Built at 2b5e0c9ec on lgv-AA-caudate-rescore-20260913 and rdg-AA-crossregion-20260918. 2654 VMRs linked to 612 loci; 84 loci with CpG-meQTL support; axis lower_in_scz_linked (significant); 145 loci with transcriptional coupling; 3 of 5 prioritized loci with GTEx support; 232 claimable ancestry-matched colocalizations. Decision 1 CAUDATE_MAGNITUDE_CLAIM_NOT_SUPPORTED (Module 08 tier 3: donor count is a plausible major contributor). Caudate is batch-confounded. |
| scz-AA-dlpfc-20260918 | AA | dlpfc | vmrset-AA-dlpfc-856067dfe289 | 2026-09-19 | Kynon J. Benjamin | PASS_SCZ_APPLICATION_QC | Built at 2b5e0c9ec on lgv-AA-dlpfc-rescore-20260913. 2063 VMRs linked to 612 loci; 58 loci with CpG-meQTL support; axis lower_in_scz_linked (significant); 99 loci with transcriptional coupling; 4 of 5 prioritized loci with GTEx support; 175 claimable ancestry-matched colocalizations. Decision 1 not applicable. |
| scz-AA-hippocampus-20260918 | AA | hippocampus | vmrset-AA-hippocampus-2d907b892215 | 2026-09-19 | Kynon J. Benjamin | PASS_SCZ_APPLICATION_QC | Built at 2b5e0c9ec on lgv-AA-hippocampus-rescore-20260913. 2017 VMRs linked to 612 loci; 64 loci with CpG-meQTL support; axis lower_in_scz_linked (significant); 66 loci with transcriptional coupling; 5 of 5 prioritized loci with GTEx support; 174 claimable ancestry-matched colocalizations. Decision 1 not applicable. |

## Contract

This module follows AGENTS.md §5.2: `_h/` holds code, `_m/` holds generated
output under immutable `runs/{RUN_ID}/` directories, `tests/` holds gitignored
smoke checks. Configuration lives in `config/` at the repository root.
