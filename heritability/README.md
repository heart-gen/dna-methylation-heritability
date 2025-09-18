# Elastic-net SNP Heritability Testing

This project aims to quantify environmental and genetic
contributions to DNA methylation in an African genetic ancestry (AA) cohort. 
To this end, we developed a comprehensive pipeline to prepare publicly
available whole genome bisulfite sequencing (WGBS), genotype, and ancestry 
data and subsequently apply a GPU-based elastic net regression software 
(https://github.com/heart-gen/GENBoostGPU) to estimate SNP heritability 
in the caudate nucleus, dorsolateral prefrontal cortex (DLPFC), and the 
hippocampus.

---

## Pipeline Overview

- [01.get_cpg_stats](#Script-`01.get_cpg_stats.R`)
- [02.pca.R](#Script-`02.pca.R`)
- [02b.res_var.R](#Script-`02b.res_var.R`)
- [03.extract_vmr.R](#Script-`03.extract_vmr.R`)
- [04.cal_vmr.sh](#Script-`04.cal_vmr.sh`)
- [05.extract_snps.sh](#Script-`05.extract_snps.sh`)

---

### Script `01.get_cpg_stats.R`

**Overview:**
- Load and filter the BSobj for selected samples
- Exclude low coverage ENCFF356LFX exclusion list regions
- Calculate methylation values for all CpG sites
- Calculate SD and mean of DNA methylation
- Write covariate files

**Inputs:**
- `BSobj.rda`
- `phenotypes-AA.tsv`
- `TOPMed_LIBD.AA.psam`

**Generated Output Files:**
- `cpg/`
    - `chr_*/`
        - `stats.rda:` SD and mean
        - `cpg_meth.phen:` Methylation values (0-1)
        - `cpg_pos.txt:` Genomic position of CpG sites
- `covs/`
    - `chr_*/`
        - `TOPMed_LIBD.AA.covar`
        - `TOPMed_LIBD.AA.qcovar`
 
 **Submission:**
 `sbatch ../_h/step_1.sh`

---

### Script `02.pca.R`

**Overview:**
- Regress out african genetic ancestry
- Perform principal component analysis on top 1 million variable CpG
- Plot principal components
- Correlate principal components with genetic ancestry

**Inputs:**
- `cpg/`
    - `chr_*/`
        - `stats.rda`
- structure.out_ancestry_proportion_raceDemo_compare

**Generated Output Files:**
- `pca/`
    - `chr_*/`
        - `pc.csv:`
        - `pca.pdf:` 
        - `pc_ances_cor.csv:`
 
 **Submission:**
 `sbatch ../_h/step_2.sh`

---

### Script `02b.res_var.R`

**Overview:**
- Chunks data for memory optimization
- Filters out CpG sites at C/T SNP sites
- Regresses out top 5 principal components generated from `02.pca.R`
- Calculates variance of residuals 

**Inputs:**
- `cpg/`
    - `chr_*/`
        - `cpg_meth.phen`
        - `cpg_pos.txt`
- `pca/`
    - `chr_*/`
        - `pc.csv:`
- `phenotypes-AA.tsv`
- `structure.out_ancestry_proportion_raceDemo_compare`
- `snps_CT`
    - `chr_*/`

**Generated Output Files:**
- `cpg/`
    - `chr_*/`
        - `res_cpg_meth.phen`
- `pca/`
    - `chr_*/`
        - `res_var_all.tsv`
 
 **Submission:**
 `sbatch ../_h/step_2b.sh`

---

### Script `03.extract_vmr.R`

**Overview:**
- Identify VMRs using residuals outputted from `02b.res_var.R`
- Write to .bed file 

**Inputs:**
- `pca/`
    - `chr_*/`
        - `res_var_all.tsv`

**Generated Output Files:**
- `vmr/`
    - `chr_*/`
        - `vmr.bed`

 **Submission:**
 `sbatch ../_h/step_3.sh`

---

### Script `04.cal_vmr.R`

**Overview:**
- Calculate methylation values for extracted VMRs
- Write to .phen file 

**Inputs:**
- `cpg/`
    - `chr_*/`
        - `stats.rda`
- `TOPMed_LIBD.AA.psam`
- `vmr.bed`: Tab-delimited file with 3 columns (chr, start, end)
             detailing all regions of interest
  
Ex.
|chr|start                        |end   |
|---|-----------------------------|------|
|1  |100565718                    |100569135|
|1  |100944974                    |100948292|
|1  |101308406                    |101308887|
|1  |104411888                    |104413497|

Generated using `cat ./vmr/chr_*/vmr.bed | sort -k1,1n > vmr.bed`

**Generated Output Files:**
- `vmr/`
    - `chr_*/`
        - `${START}_${END}_meth.phen`

 **Submission:**
 `sbatch ../_h/step_4.sh`

---

### Script `05.extract_snps.sh`

**Overview:**
- Use Plink to extract SNPs 500 kb around each VMR
- Extract subset for samples used in study
- Check for potential window selection errors 

**Inputs:**
- `TOPMed_LIBD.AA.pgen`
- `TOPMed_LIBD.AA.psam`
- `TOPMed_LIBD.AA.pvar`
- `chromosome_sizes.txt`
- `vmr.bed`
- `samples.txt`

**Generated Output Files:**
- `plink_format/`
    - `chr_*/`
        - `TOPMed_LIBD.AA.${START}_${END}`
 
 **Submission:**
 `sbatch ../_h/step_5.sh`
 
