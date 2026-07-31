# Public brain meQTL resources: source and download notes

This document accompanies `inputs/data_dictionary/_m/public_meqtl_resources.tsv`.

## Inclusion criteria (from analysis_strategy.md)

Prefer resources with:

- human postmortem brain tissue
- CpG-level methylation measurements
- cis-meQTL summary results or lead meQTL lists
- known (or ascertainable) genome build
- clear brain region
- adequate sample size
- accessible ancestry information
- sufficient overlap with WGBS CpGs in VMRs

Do not pool incompatible platforms without justification.

## Primary resource

### brainseq_wgbs_meqtl

- **Citation:** Perzel Mandell AM et al. Genome-wide sequencing-based identification of methylation quantitative trait loci and their role in schizophrenia risk. *Nature Communications* 2021.
- **Why primary:** WGBS in adult DLPFC and hippocampus — platform-matched to this study.
- **Browser:** https://eqtl.brainseq.org/WGBS_meQTL/
- **Download plan:** full catalogs on Synapse `syn25992404` (DOI 10.7303/syn25992404) under PsychENCODE controlled access — **not available for this project currently**. Nature supplements only contain SCZ-risk SNP–CpG tables (interim `brainseq_wgbs_meqtl_scz_subset`).
- **Harmonization:** WGBS paper used GRCh38; do not require exact lead-SNP identity for enrichment tests.
- **Decision:** treat BrainSeq published full catalog as deferred. Phase 3 external validation proceeds with Jaffe (DLPFC 450K) and Schulz (hippocampus array). Regenerating meQTLs on overlapping BrainSeq AA/EA discovery samples is **not** independent external validation; that role is already filled by Phase 1 internal cis-meQTL mapping.

## Secondary resources

### jaffe_dlpfc_450k_meqtl

- Adult DLPFC Illumina 450K cis-meQTL catalogs from LIBD / Jaffe et al.
- Useful positive-control overlap for array-detectable CpGs within VMRs.
- Expect moderate CpG overlap with WGBS; analyze separately from WGBS external tests.

### schulz_hippocampus_array_meqtl

- Hippocampal array cis-meQTL catalog (Nat Commun 2017).
- Supplementary data / sciebo link in the paper.
- Region-matched external support for hippocampus VMRs.

## Exploratory only

### hannon_fetal_brain_meqtl

- Fetal brain context; not adult postmortem.
- Optional positive-control overlap only; not a primary external validation resource.

## Download checklist (Phase 3)

1. Create `meqtl-validation/03_external_meqtl_validation/_m/raw/{resource_id}/`
2. Save original files + checksums / URLs / download date in `download_manifest.tsv`
3. Run harmonization scripts to produce `*_hg38_harmonized.tsv.gz`
4. Quantify CpG overlap with VMR CpGs per region
5. Run external support ~ predictability models using `config/analysis_thresholds.yml`
