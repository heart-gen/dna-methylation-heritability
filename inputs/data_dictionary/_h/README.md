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
- **Download plan:** retrieve cis-meQTL summary tables from BrainSeq / PsychENCODE release associated with the paper; store under `meqtl-validation/03_external_meqtl_validation/_m/raw/brainseq_wgbs/` without overwriting source archives.
- **Harmonization:** confirm genome build on ingest; map CpG coordinates to hg38; record liftOver success/failure if needed; do not require exact lead-SNP identity for enrichment tests.

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
