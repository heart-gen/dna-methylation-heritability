# GREML Heritability Testing

This project aims to quantify environmental and genetic
contributions to DNA methylation in an African genetic ancestry (AA) cohort. 
To this end, we developed a comprehensive pipeline to prepare publicly
available whole genome bisulfite sequencing (WGBS), genotype, and ancestry data 
and subsequently apply a genomic restricted maximum likelihood (GREML) mixed model 
to estimate SNP heritability in the hippocampus nucleus. This pipeline leverages 
the Genome-wide Complex Trait Ananlyis [`(GCTA)`](https://github.com/jianyangqt/gcta) 
software for heritability testing.  

---

## Pipeline Overview

- [01.extract_vmr.R](#Script-`01.extract_vmr.R`)
- [02.pca.R](#Script-`02.pca.R`)
- [03.cal_vmr.R](#Script-`03.cal_vmr.R`)
- [04.extract_snps.sh](#Script-`04.extract_snps.sh`)
- [05.greml.sh](#Script-`05.greml.sh`)
- [06.summary.py](#Script-`06.summary.py`)

---

### Script `01.extract_vmr.R`

**Overview:**
- Load and filter the BSobj for selected sample
- Calculate methylation values for all CpG sites
- Calculate SD and mean of DNA methylation
- Extract variably methylated regions
- Write covariate files 

**Inputs:**
- `hippocampus_chr*_BSobj.rda`
- `merged_phenotypes.csv`
- `TOPMed_LIBD.AA.psam`

**Generated Output Files:**
- `cpg/`
    - `chr_*/`
        - `stats.rda:` SD and mean
        - `cpg_meth.phen:` Methylation values (0-1) 
- `vmr/`
    - `chr_*/`
        - `vmr.bed`
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

### Script `03.cal_vmr.R`

**Overview:**
- Calculate methylation values for extracted VMR
- Write to .phen file 

**Inputs:**
- `cpg/`
    - `chr_*/`
        - `stats.rda`
- `TOPMed_LIBD.AA.psam`
- `vmr_list.txt`: Tab-delimited file with 3 columns (chr, start, end)
                  detailing all regions of interest
  
Ex.
|chr|start                        |end   |
|---|-----------------------------|------|
|1  |100565718                    |100569135|
|1  |100944974                    |100948292|
|1  |101308406                    |101308887|
|1  |104411888                    |104413497|

Generated using `cat ./chr_*/vmr.bed | sort -k1,1n > vmr_list.txt`

**Generated Output Files:**
- `vmr/`
    - `chr_*/`
        - `${START}_${END}_meth.phen`

 **Submission:**
 `sbatch ../_h/step_3.sh`

---

### Script `04.extract_snps.sh`

**Overview:**
- Use Plink to extract SNPs 500 kb around each VMR
- Check for potential window selection errors 

**Inputs:**
- `TOPMed_LIBD.AA.pgen`
- `TOPMed_LIBD.AA.psam`
- `TOPMed_LIBD.AA.pvar`
- `chromosome_sizes.txt`
- `vmr_list.txt`

**Generated Output Files:**
- `plink_format/`
    - `chr_*/`
        - `TOPMed_LIBD.AA.${START}_${END}`
 
 **Submission:**
 `sbatch ../_h/step_4.sh`

---

### Script `05.greml.sh`

**Overview:**
- Calculate SNP LD scores and stratify based on individual SNP scores
- Plot histogram of SNP LD scores for each region
- Create genetic relation matricies (GRM) for each stratified group 
- Perform GREML-LDMS analysis using multiple GRM

**Inputs:**
- `../_h/05.stratify_LD.R`
- `plink_format/`
    - `chr_*/`
        - `TOPMed_LIBD.AA.${START}_${END}`
- `covs/`
    - `chr_*/`
        - `TOPMed_LIBD.AA.covar`
        - `TOPMed_LIBD.AA.qcovar`
- `vmr/`
    - `chr_*/`
        - `${START}_${END}_meth.phen`
- `vmr_list.txt`

**Generated Output Files:**
- `h2/`
    - `chr_*/`
        - `TOPMed_LIBD.AA.${START}_${END}.score.ld`
        - `TOPMed_LIBD.AA.${START}_${END}_group{1-4}.grm.bin`
        - `TOPMed_LIBD.AA.${START}_${END}_group{1-4}.grm.id`
        - `TOPMed_LIBD.AA.${START}_${END}_group{1-4}.grm.N.bin`
        - `TOPMed_LIBD.AA.${START}_${END}_group{1-4}.log`
        - `TOPMed_LIBD.AA.${START}_${END}.log`
        - `TOPMed_LIBD.AA.${START}_${END}.hsq`
        - `${START}_${END}_multi_GRMs.txt`
        - `hist/`
            - `${START}_${END}_ld_hist.pdf`
          
 **Submission:**
 `sbatch ../_h/step_5.sh`

---

### Script `06.summary.py`

**Overview:**
- Write GCTA output to .csv files
- Compile results for all regions

**Inputs:**
- `vmr_list.txt`
- `h2/`
    - `chr_*/`
        - `TOPMed_LIBD.AA.${START}_${END}.hsq`

**Generated Output Files:**
- `h2/`
    - `chr_*/`
        - `TOPMed_LIBD.AA.${START}_${END}.csv`
- `summary/`
    - `TOPMed_LIBD.AA.csv`

 **Submission:**
 `sbatch ../_h/step_6.sh`
 
