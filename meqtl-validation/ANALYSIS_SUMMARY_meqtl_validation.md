## CpG cis-meQTL evidence supports genetically anchored methylation variability in human brain

### Purpose

Variably methylated regions (VMRs) identify genomic intervals where DNA
methylation differs substantially among individuals. The manuscript uses an
elastic-net score to rank VMRs by how well nearby common SNPs predict their
methylation. Biologically, a high score should mark a VMR whose interindividual
methylation variability is anchored by local genetic regulation. Prediction
alone, however, does not show that individual CpGs in that VMR have conventional
cis-meQTLs, and selected elastic-net SNPs are not themselves significant or
causal meQTLs.

This analysis therefore asked whether VMRs with greater local genetic
predictability contain a larger proportion of CpGs with independently mapped
cis-meQTL evidence. It also asked whether this relationship is observed across
caudate, dorsolateral prefrontal cortex (DLPFC), and hippocampus; whether it is
present in both Black American and non-Hispanic white American donor groups;
whether meQTL-supported VMRs are more often coupled to expression or splicing;
and whether genetically anchored VMRs remain concentrated in repeat-rich and
repressive genomic compartments after technical and cell-composition
sensitivity analyses.

The biological interpretation is deliberately focused: a positive association
between local genetic predictability and CpG cis-meQTL burden would validate the
VMR ranking as a measure of genetically anchored methylation architecture. It
would not establish locus-level heritability, causal variants, methylation
mediation, or ancestry-specific effects.

### Current analysis scope

This summary reports only outputs regenerated from the final CpG and VMR tables
using analysis schema version 2. The current evidence set includes internal
CpG-level mapping, VMR-level burden models, donor-group analyses,
expression/splicing integration, repeat and mappability sensitivity, cell-type
sensitivity, and the external-resource audit. Cross-region sharing,
method-matched caudate downsampling, the schizophrenia-risk application, and the
simulation-calibrated estimator sensitivity have not yet been rerun on the same
final post-QC inputs and are not interpreted here.

| Analysis component | Current status | Manuscript use |
|---|---|---|
| Internal Black American CpG cis-meQTL mapping | Final post-QC outputs available | Primary evidence |
| VMR predictability–meQTL burden | Final schema-v2 outputs available | Primary evidence |
| White American stratified mapping and donor-group comparison | Final schema-v2 outputs available | Portability analysis |
| Expression and splicing integration | Regenerated from final burden tables | Biological interpretation |
| Repeat, mappability, and cell-type sensitivity | Regenerated from final burden tables | Technical robustness |
| Independent public brain meQTL gradient | Audited, but not estimable from available tables | Unmet external criterion |
| Cross-region sharing and caudate downsampling | Final-input rerun required | Do not cite yet |
| Schizophrenia-risk application | Final-input rerun required | Do not cite yet |
| Calibrated local SNP-variance sensitivity | Final-input rerun required | Do not cite yet |

### Inputs

- Predefined hg38 VMR intervals and their elastic-net local genetic-predictability
  scores.
- Raw, non-residualized WGBS CpG methylation values and imputed genotype dosages
  from matched BrainSeq donors.
- Black American discovery samples: caudate N = 153, DLPFC N = 111, and
  hippocampus N = 116.
- Non-Hispanic white American stratified samples: caudate N = 129, DLPFC N = 55,
  and hippocampus N = 60.
- The locked M3a discovery covariates: age at death, sex, schizophrenia
  diagnosis, genotype principal components 1–5, and methylation principal
  components 1–5.
- Existing VMR–expression and VMR–splicing association results; these analyses
  were reused rather than initiating a new transcriptome-wide screen.
- LINE/L1, H3K9me3, quiescent-chromatin, mappability, segmental-duplication,
  SNP-proximity, CpG, coverage, and broad genomic annotations.
- RNA MuSiC cell-composition estimates in all three regions and DNAm scMD
  estimates for the prespecified caudate sensitivity.
- Jaffe DLPFC 450K [@doi:10.1038/nn.4181], Schulz hippocampal array
  [@doi:10.1038/s41467-017-01818-4], and BrainSeq WGBS
  [@doi:10.1038/s41467-021-25517-3] meQTL resources.

### Methods Text

We mapped CpG-level cis-meQTLs separately in caudate, DLPFC, and hippocampus.
The primary molecular phenotype was the non-residualized methylation level of an
individual CpG contained within a predefined VMR. CpGs were required to have
coverage of at least five reads in at least 80% of samples and nonzero variance.
Blacklisted CpGs and CpGs directly overlapping common C/T SNPs were excluded
before phenotype construction. Genotypes were restricted to autosomal variants
with minor-allele frequency at least 0.05, imputation quality at least 0.8,
missingness at most 0.05, and Hardy–Weinberg P at least 1 × 10^-6.

TensorQTL tested SNPs within ±500 kb of each CpG using the locked M3a covariate
model. Storey q-values were calculated separately within each brain region, and
a CpG was classified as meQTL-supported when its lead cis association had
q ≤ 0.05. The non-Hispanic white American stratified analysis used the available
baseline demographic, diagnostic, and genotype-PC covariates and was treated as
a secondary portability analysis rather than a direct power-matched replication
of M3a.

For each coordinate-resolved VMR, CpG meQTL burden was defined as the number of
meQTL-supported CpGs divided by the number of CpGs with an actual TensorQTL test
result. CpGs prepared for mapping but lacking a test result were retained in QC
counts but were not coded as negative evidence. The primary VMR model was an
overdispersion-scaled binomial generalized linear model with HC3 robust
covariance, fitted separately by region. The standardized continuous local
genetic-predictability score was the primary predictor. The prespecified model
adjusted, where sufficiently complete, for tested CpG number, mean coverage,
methylation mean and variance, VMR length, CpG density, local tested-SNP
opportunity, mappability, LINE/L1 overlap, problematic-region overlap, and broad
genomic annotation.

As a complementary extreme-group analysis, VMRs in the highest and lowest
predictability quintiles were propensity-score matched without replacement
using a 0.25-SD caliper and exact matching on broad genomic annotation.
Within-pair label randomization with 10,000 permutations tested the difference
in meQTL burden. A matched result was considered interpretable only when the
maximum absolute post-match standardized mean difference was no greater than
0.10.

Donor-group portability was assessed using the correlation of VMR
predictability scores and by comparing identical lead SNP–CpG pairs across
stratified analyses. MAF and cis-SNP testing opportunity were matched before
interpreting discovery-rate differences. Expression and splicing analyses
modeled whether a VMR with any CpG meQTL support was more likely to have an
existing significant VMR–transcript association, adjusting for the number of
tested features, VMR length, VMR-to-feature distance, methylation variance, and
local SNP number.

Genomic-compartment sensitivity analyses modeled LINE/L1, H3K9me3, and
quiescent-chromatin overlap as functions of standardized VMR predictability and
repeated the analyses in high-mappability intervals and after excluding
SNP-proximal or segmental-duplication intervals. Cell-type sensitivity adjusted
for the proportion of each VMR's methylation variance associated with RNA MuSiC
cell-composition PCs. Caudate DNAm scMD PCs were included as a prespecified
orthogonal sensitivity; DLPFC and hippocampal DNAm estimates remained
exploratory because they did not pass the DNAm–RNA integration gate.

External enrichment was considered estimable only when a public resource
provided both supported and unsupported CpGs or VMRs within a documented assayed
universe. VMRs absent from a published positive list were not treated as
external negatives. BrainSeq WGBS results were not eligible as independent
replication because the donor cohort overlaps the present study.

### Results Text

**CpG-level mapping identified extensive local genetic regulation of brain
methylation.** In the Black American discovery samples, 195,830 caudate CpGs,
152,276 DLPFC CpGs, and 154,600 hippocampal CpGs had a valid cis test. At
region-specific q ≤ 0.05, 80,281 caudate CpGs (41.0%), 50,288 DLPFC CpGs
(33.0%), and 54,898 hippocampal CpGs (35.5%) had cis-meQTL support. These counts
show that nearby common variation contributes to methylation differences at a
large fraction of variable brain CpGs.

Genome-wide lambda values calculated from the lead-cis P values were 4.96,
3.71, and 3.91 in caudate, DLPFC, and hippocampus. Because this statistic mixes
true dense cis signal with residual calibration, it should not be interpreted as
a null-only inflation estimate. Among CpGs not passing the regional FDR
threshold, lambda was lower but remained above one (1.45, 1.46, and 1.41), so
the large discovery counts should be presented together with the locked
covariate-model QC rather than as evidence of perfect null calibration.

**Greater VMR local genetic predictability was strongly associated with greater
CpG cis-meQTL burden in every brain region.** The schema-v2 aggregation included
195,830 tested CpGs across 11,372 caudate VMRs, 152,276 CpGs across 9,975 DLPFC
VMRs, and 154,600 CpGs across 9,800 hippocampal VMRs. Mean VMR meQTL burden was
0.42, 0.31, and 0.34, respectively. The complete-case prespecified models
retained 8,738 caudate, 3,962 DLPFC, and 7,152 hippocampal VMRs. After technical
and genomic adjustment, each standard-deviation increase in local predictability was
associated with 22.57-fold greater odds that a tested CpG had meQTL support in
caudate (95% CI 17.73–28.74; P = 3.3 × 10^-141), 6.50-fold greater odds in DLPFC
(95% CI 4.85–8.71; P = 5.5 × 10^-36), and 8.60-fold greater odds in hippocampus
(95% CI 7.19–10.28; P = 1.4 × 10^-123). Thus, the manuscript's VMR ranking is
strongly aligned with a conventional CpG-level measure of local genetic
regulation and is not simply detecting VMR size, CpG density, coverage, local
SNP opportunity, or the measured genomic annotations.

The extreme-group contrasts pointed in the same biological direction. Matched
high-predictability VMRs had mean meQTL-burden differences of 0.73 across 679
caudate pairs, 0.53 across 251 DLPFC pairs, and 0.63 across 523 hippocampal pairs
(all permutation P = 1.0 × 10^-4).
Post-match balance passed in caudate (maximum absolute SMD = 0.087) but not in
DLPFC (0.179) or hippocampus (0.192). The adjusted continuous models therefore
provide the primary three-region evidence. The matched result is confirmatory
in caudate and directionally supportive, but not balance-qualified, in the other
two regions.

**The aggregate genetic-predictability gradient was also present in the white
American stratified analyses.** Of 196,740, 172,980, and 148,968 tested CpGs,
70,678 caudate, 12,738 DLPFC, and 31,171 hippocampal CpGs had meQTL support.
Adjusted predictability coefficients were positive in all three regions
(log-odds coefficients = 3.43, 1.12, and 1.38; all P ≤ 1.4 × 10^-89). VMR
predictability scores were most portable in caudate (Spearman rho = 0.51),
moderately portable in hippocampus (rho = 0.30), and weakly correlated in DLPFC
(rho = 0.12). Among comparisons restricted to the identical lead variant–CpG
pair, direction concordance for pairs significant in both donor groups was
0.836, 0.824, and 0.816.

Matching on MAF and cis-SNP testing opportunity eliminated the discovery-rate
difference in caudate and DLPFC but not hippocampus. Together, these results
support portability of the aggregate genetically anchored architecture while
showing that locus-level detectability varies with allele frequency, testing
opportunity, sample size, LD, or effect heterogeneity. They do not support a
claim of ancestry-specific methylation biology.

**CpG meQTL-supported VMRs were more likely to be transcriptionally coupled.**
After adjustment for VMR and feature-testing characteristics, meQTL-supported
VMRs had greater odds of an existing expression association in caudate
(OR = 1.81, P = 1.5 × 10^-6), DLPFC (OR = 12.21,
P = 3.2 × 10^-11), and hippocampus (OR = 8.76, P = 8.9 × 10^-9). Splicing
associations were enriched in caudate (OR = 4.55, P = 9.5 × 10^-15) and DLPFC
(OR = 1.77, P = 0.013), but not hippocampus (OR = 0.72, P = 0.55). These results
place genetically anchored VMRs in a broader regulatory context: VMRs with
direct CpG cis-meQTL evidence are more likely to show detectable relationships
with transcription, especially gene expression. The associations do not
establish that methylation mediates genotype effects on expression or splicing.

**Genetically anchored VMRs were concentrated in repressive chromatin, with a
regionally qualified LINE/L1 result.** In adjusted models, greater local
predictability was associated with H3K9me3 overlap in caudate, DLPFC, and
hippocampus (OR = 1.77, 1.32, and 1.52) and with quiescent chromatin in all three
regions (OR = 2.47, 1.38, and 1.72; all P ≤ 3.0 × 10^-9). Direction was retained
after high-mappability, SNP-proximity, and segmental-duplication sensitivity
analyses. LINE/L1 overlap was also enriched in the adjusted models
(OR = 1.84, 1.19, and 1.32), but the DLPFC estimate reversed and became
nonsignificant after the high-mappability restriction. The defensible
interpretation is therefore that LINE/L1 enrichment is robust in at least two
regions, not uniformly brain-wide.

Adjustment for cell-composition–correlated methylation properties preserved the
predictability–meQTL burden association in all three regions. H3K9me3 and
quiescent-chromatin associations also remained positive and significant in all
regions. LINE/L1 enrichment remained significant in caudate (OR = 1.56) and
hippocampus (OR = 1.19) but not DLPFC (OR = 1.05). The caudate findings were also
directionally retained with DNAm-derived cell PCs. These results reduce concern
that the repeat and repressive-compartment patterns are explained entirely by
bulk cellular mixture, while recognizing that cell-type adjustment cannot
recover cell-specific meQTL effects from bulk tissue.

**The available public tables document overlap with published brain meQTLs but
cannot test the required external gradient.** All 3,815 Jaffe DLPFC CpGs
represented in the current DLPFC VMR overlap table were labeled supported,
leaving no assayed unsupported comparison group. The Schulz hippocampal table
was likewise positive-only for all 966 overlapping hippocampal CpGs.
Consequently, neither table can estimate whether the
probability of external support increases with VMR predictability. BrainSeq WGBS
results are biologically relevant but not independent because of cohort
overlap. The prespecified independent external-validation criterion is therefore
not met; this is an absence of an estimable external contrast, not evidence that
the strong internal gradient is false.

**Overall interpretation.** The current analysis provides strong internal
validation that the VMR local genetic-predictability axis identifies genetically
anchored brain methylation. The relationship is large, persists after extensive
technical adjustment in all three regions, is observed in both donor groups,
and connects to transcriptional and repressive-chromatin architecture. The
strict full-validation gate is not yet satisfied because the public-resource
gradient is not estimable and the prespecified matched-balance criterion passed
only in caudate. The manuscript can state that CpG cis-meQTL burden independently
supports the VMR predictability gradient, while describing external validation
and two-region matching as unresolved rather than complete.

### Figure and Table Notes

- Potential main figure: internal CpG cis-meQTL validation of the VMR
  predictability gradient.
  - Panel A: region-specific distributions of VMR meQTL burden across
    predictability quantiles.
  - Panel B: adjusted odds ratios per standard-deviation increase in continuous
    predictability, with 95% confidence intervals.
  - Panel C: high-versus-low matched burden differences, displaying post-match
    maximum absolute SMD and marking only caudate as balance-qualified.
  - Key message: more predictable VMRs contain substantially more CpGs with
    conventional cis-meQTL support in all three regions.

- Potential second main or supplementary figure: biological context of
  genetically anchored VMRs.
  - Panel A: adjusted expression and splicing enrichment for meQTL-supported
    VMRs.
  - Panel B: LINE/L1, H3K9me3, and quiescent-chromatin estimates across technical
    restrictions.
  - Panel C: MuSiC cell-composition sensitivity, with caudate DNAm scMD as an
    orthogonal sensitivity.
  - Key message: genetically anchored VMRs are transcriptionally coupled and
    enriched in repressive compartments; LINE/L1 enrichment is regionally
    heterogeneous.

- Potential supplementary tables:
  - CpG mapping and calibration QC by region and donor group.
  - Complete adjusted VMR burden models and matched-balance diagnostics.
  - Identical-pair donor-group concordance and MAF/testing-opportunity matching.
  - Expression/splicing enrichment models.
  - Consolidated repeat, mappability, problematic-region, and cell-type
    sensitivity estimates.
  - External-resource assay-universe audit showing why the gradient is not
    estimable.

A consolidated manuscript figure has not yet been regenerated from the final
post-QC inputs. Pre-final cross-region, downsampling, calibrated-estimator, and
schizophrenia locus figures should not be used for manuscript claims.

### Reproducibility Information

- Analysis root: `meqtl-validation/`.
- Primary configuration: `config/meqtl_parameters.yml`,
  `config/covariates.yml`, `config/analysis_thresholds.yml`, and
  `config/paths.yml`.
- Primary CpG-mapping scripts:
  `01_cpg_meqtl_mapping/_h/02a_prepare_cpg_bed.py`,
  `02b_merge_cpg_beds.py`,
  `03_prepare_genotypes.sh`, `04_tensorqtl_map.py`, and
  `05_qc_summarize.py`.
- VMR aggregation and inference:
  `02_vmr_meqtl_burden/_h/01_aggregate_vmr_burden.py` and
  `02_fit_burden_models.py`.
- Donor-group, transcription, repeat, and cell-type modules:
  `05_donor_group_comparison/_h/`,
  `06_transcription_splicing_integration/_h/`,
  `07_repeat_mappability_sensitivity/_h/`, and
  `11_celltype_compartment_sensitivity/_h/`.
- SLURM execution: TensorQTL mapping uses the `gengpu` partition with one A100
  GPU; downstream models use the `genomics` partition. Each module has
  repository-local `step_*.sh` submission scripts. There is not yet a single
  end-to-end submission wrapper for the complete general meQTL-validation
  workflow.
- Compute environment: conda environment
  `/projects/p32505/opt/envs/genomics`; no container was used.
- Environment at summary time: Python 3.11.13, TensorQTL 1.0.10, NumPy 2.2.6,
  and pandas 2.3.3. A run-frozen project-local environment specification is not
  currently attached to the general workflow.
- Primary random seed: `20260722`; matched analyses used 10,000 within-pair
  permutations.
- Current output refresh date: 2026-08-08.
- Repository commit at summary time: `e94c103b9` with additional uncommitted
  analysis outputs and documentation changes.
- Analysis schema: version 2; downstream consumers reject older burden tables.

Key current outputs are:

- `01_cpg_meqtl_mapping/{region}/_m/tensorqtl/qc/meqtl_qc_summary.tsv`
- `01_cpg_meqtl_mapping/{region}/_m/tensorqtl/qc/calibration/lambda_by_qval.tsv`
- `02_vmr_meqtl_burden/_m/{region}/aggregation_summary.tsv`
- `02_vmr_meqtl_burden/_m/{region}/burden_model_results.tsv`
- `02_vmr_meqtl_burden/_m/{region}/matched_analysis_results.tsv`
- `05_donor_group_comparison/_m/aa_ea_predictability_portability.tsv`
- `05_donor_group_comparison/_m/aa_ea_effect_concordance_summary.tsv`
- `05_donor_group_comparison/_m/maf_ld_matched_discovery_claim.tsv`
- `06_transcription_splicing_integration/_m/tx_enrichment_primary.tsv`
- `07_repeat_mappability_sensitivity/_m/consolidated_robustness_table.tsv`
- `11_celltype_compartment_sensitivity/_m/celltype_claim_summary.tsv`
- `03_external_meqtl_validation/_m/phase3_criterion5_verdict.tsv`

### Limitations and Integration Notes

The established elastic-net score should continue to be called local genetic
predictability or aggregate local SNP contribution. Although it is strongly
validated by CpG-level meQTL burden, it is not a calibrated locus-level
heritability estimate. Elastic-net-selected SNPs are neither significant meQTLs
nor fine-mapped causal variants.

The large lead-cis lambda values reflect a mixture of abundant true local signal
and residual calibration; they should not be used alone to argue that all test
statistics are well calibrated. The M3a choice and nonsignificant-site lambda
values should remain visible in supplementary QC.

The adjusted continuous result is the primary validation because DLPFC and
hippocampal extreme-group matches did not meet the prespecified balance
threshold. Technical-covariate completeness also reduced the adjusted DLPFC
analysis to 3,962 VMRs, so its effect estimate applies to that complete-case
subset. The donor-group analysis is stratified but not power matched, and the
white American DLPFC sample is small. Discovery differences must therefore not
be interpreted as ancestry-specific without adequately powered interaction
tests and LD-aware comparisons.

Expression and splicing enrichment establishes regulatory coupling, not
mediation or causality. The hippocampal splicing result is null and should remain
visible. LINE/L1 enrichment is sensitive to the DLPFC high-mappability
restriction; H3K9me3 and quiescent-chromatin results are more consistent across
regions. Bulk cell-composition sensitivity reduces a major confounding concern
but does not identify the cell type in which a cis-meQTL acts.

Independent external comparative validation remains the principal evidence gap.
The current public tables cannot supply a valid supported-versus-unsupported
assayed universe, and overlapping BrainSeq results cannot fill that role. A new
external dataset must provide genome-wide or array-wide tested CpGs, including
nonsignificant results, before the external gradient can be estimated.

Before final manuscript packaging, rerun the cross-region sharing,
rate-based caudate downsampling, schizophrenia-risk application, and calibrated
estimator sensitivity from the current final CpG and schema-v2 VMR inputs. Then
regenerate the consolidated figure, supplementary tables, and the Phase 1 lock
memo without retaining superseded result versions.

### Claim Checklist

| Claim | Current support |
|---|---|
| Higher VMR local predictability is associated with greater CpG cis-meQTL burden | Supported in 3/3 regions after prespecified adjustment |
| The association is independent of measured technical and genomic features | Supported in 3/3 regions |
| The extreme-group matched result is fully balanced | Supported in caudate only; DLPFC and hippocampus are directional only |
| The aggregate gradient is present in both donor groups | Supported in 3/3 regions |
| Donor-group differences are ancestry-specific | Not supported and not claimed |
| meQTL-supported VMRs are enriched for expression associations | Supported in 3/3 regions |
| meQTL-supported VMRs are enriched for splicing associations | Supported in caudate and DLPFC; null in hippocampus |
| H3K9me3 and quiescent enrichment survive technical and cell-type sensitivity | Supported in 3/3 regions |
| LINE/L1 enrichment is uniformly robust across regions | Not supported; robust in at least two regions, DLPFC is fragile |
| An independent public brain dataset validates the predictability gradient | Not currently estimable |
| Cross-region, caudate-downsampling, schizophrenia, and calibrated-estimator claims are final | No; final-input reruns are required |
