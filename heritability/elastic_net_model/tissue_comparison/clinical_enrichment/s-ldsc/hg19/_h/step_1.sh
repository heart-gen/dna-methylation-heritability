#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --time=02:00:00
#SBATCH --ntasks-per-node=8
#SBATCH --job-name=munge_sumstats
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --output=logs/munge_stats.%j.log

# =============================================================================
# Step 1: Munge GWAS Summary Statistics
# =============================================================================
# This script processes GWAS summary statistics for all diseases using
# LDSC's munge_sumstats.py to prepare them for S-LDSC analysis.
# =============================================================================

# Source configuration
SCRIPT_DIR="/ocean/projects/bio250020p/kbenjamin/projects/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/s-ldsc/hg19/_h"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/genomics

# Create output directory
OUT_DIR="./sumstats"
mkdir -p "$OUT_DIR"

# LDSC wrapper (patches pandas/Python 3 compatibility issues)
LDSC_WRAPPER="${SCRIPT_DIR}/ldsc_wrapper.py"

# GWAS directory base
GWAS_BASE="/ocean/projects/bio250020p/shared/resources/gwas"

# -----------------------------------------------------------------------------
# Munge each GWAS file
# -----------------------------------------------------------------------------
# Each GWAS has different column formats, so we handle them individually

# --- Alzheimer's Disease (AD) ---
# Source: PGC ALZ2 (Kunkle et al. 2019)
# Columns: chr, PosGRCh37, testedAllele, otherAllele, z, p, N
log_message "Processing: Alzheimer's Disease (ad)"
python "$LDSC_WRAPPER" "$LDSC_DIR" munge_sumstats.py \
    --sumstats "${GWAS_BASE}/PGC/AD/data/PGCALZ2sumstatsExcluding23andMe.txt.gz" \
    --out "${OUT_DIR}/ad" \
    --a1 testedAllele \
    --a2 otherAllele \
    --signed-sumstats z,0 \
    --p p \
    --N-col N \
    --merge-alleles "${RESOURCE_DIR}/hm3_no_MHC.list.txt" \
    --chunksize 500000

# --- Schizophrenia (SCZ) ---
# Source: PGC3 (Trubetskoy et al. 2022 Nature)
# Columns: CHROM, ID, POS, A1, A2, FCAS, FCON, IMPINFO, BETA, SE, PVAL, NCAS, NCON, NEFF
log_message "Processing: Schizophrenia (scz)"
python "$LDSC_WRAPPER" "$LDSC_DIR" munge_sumstats.py \
    --sumstats "${GWAS_BASE}/PGC/SCZ/PGC3/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz" \
    --out "${OUT_DIR}/scz" \
    --a1 A1 \
    --a2 A2 \
    --signed-sumstats BETA,0 \
    --p PVAL \
    --snp ID \
    --frq FCAS \
    --N-cas-col NCAS \
    --N-con-col NCON \
    --merge-alleles "${RESOURCE_DIR}/hm3_no_MHC.list.txt" \
    --chunksize 500000

# --- Major Depressive Disorder (MDD) ---
# Source: PGC MDD3 (Giannakopoulou et al. 2021)
# Columns: MarkerName, chr, pos, Allele1, Allele2, Freq1, Effect, StdErr, P.SE
log_message "Processing: Major Depressive Disorder (mdd)"
python "$LDSC_WRAPPER" "$LDSC_DIR" munge_sumstats.py \
    --sumstats "${GWAS_BASE}/mdd/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe.txt.gz" \
    --out "${OUT_DIR}/mdd" \
    --a1 Allele1 \
    --a2 Allele2 \
    --signed-sumstats Effect,0 \
    --p P.SE \
    --snp MarkerName \
    --frq Freq1 \
    --merge-alleles "${RESOURCE_DIR}/hm3_no_MHC.list.txt" \
    --chunksize 500000

# --- Bipolar Disorder (BIP) ---
# Source: PGC BIP 2021 (Mullins et al. 2021)
# Columns: #CHROM, POS, ID, A1, A2, BETA, SE, PVAL, NGT, FCAS, FCON, IMPINFO, NEFFDIV2, NCAS, NCON
log_message "Processing: Bipolar Disorder (bip)"
python "$LDSC_WRAPPER" "$LDSC_DIR" munge_sumstats.py \
    --sumstats "${GWAS_BASE}/bip/pgc-bip2021-all.vcf.tsv.gz" \
    --out "${OUT_DIR}/bip" \
    --a1 A1 \
    --a2 A2 \
    --signed-sumstats BETA,0 \
    --p PVAL \
    --snp ID \
    --frq FCAS \
    --N-cas-col NCAS \
    --N-con-col NCON \
    --merge-alleles "${RESOURCE_DIR}/hm3_no_MHC.list.txt" \
    --chunksize 500000

# --- Parkinson's Disease (PD) ---
# Source: UK Biobank imputed (self-reported)
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Parkinson's Disease (pd)"
python "$LDSC_WRAPPER" "$LDSC_DIR" munge_sumstats.py \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_UKB_20002_1262_self_reported_parkinsons_disease.txt.gz" \
    --out "${OUT_DIR}/pd" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/hm3_no_MHC.list.txt" \
    --chunksize 500000

# --- Multiple Sclerosis (MS) ---
# Source: IMMUNOBASE
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Multiple Sclerosis (ms)"
python "$LDSC_WRAPPER" "$LDSC_DIR" munge_sumstats.py \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_IMMUNOBASE_Multiple_sclerosis_hg19.txt.gz" \
    --out "${OUT_DIR}/ms" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/hm3_no_MHC.list.txt" \
    --chunksize 500000

# --- Rheumatoid Arthritis (RA) ---
# Source: OKADA trans-ethnic
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Rheumatoid Arthritis (ra)"
python "$LDSC_WRAPPER" "$LDSC_DIR" munge_sumstats.py \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_RA_OKADA_TRANS_ETHNIC.txt.gz" \
    --out "${OUT_DIR}/ra" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/hm3_no_MHC.list.txt" \
    --chunksize 500000

# --- Asthma ---
# Source: GABRIEL consortium
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Asthma (asthma)"
python "$LDSC_WRAPPER" "$LDSC_DIR" munge_sumstats.py \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_GABRIEL_Asthma.txt.gz" \
    --out "${OUT_DIR}/asthma" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/hm3_no_MHC.list.txt" \
    --chunksize 500000

# --- Coronary Artery Disease (CAD) ---
# Source: CARDIoGRAM C4D
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Coronary Artery Disease (cad)"
python "$LDSC_WRAPPER" "$LDSC_DIR" munge_sumstats.py \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.txt.gz" \
    --out "${OUT_DIR}/cad" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/hm3_no_MHC.list.txt" \
    --chunksize 500000

# --- Hypertension (HTN) ---
# Source: ICBP Systolic Blood Pressure
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Hypertension (htn)"
python "$LDSC_WRAPPER" "$LDSC_DIR" munge_sumstats.py \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_ICBP_SystolicPressure.txt.gz" \
    --out "${OUT_DIR}/htn" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/hm3_no_MHC.list.txt" \
    --chunksize 500000

# --- Stroke ---
# Source: ISGC METASTROKE (Malik et al. 2016) - All strokes
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Stroke (stroke)"
python "$LDSC_WRAPPER" "$LDSC_DIR" munge_sumstats.py \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_ISGC_Malik_2016_METASTROKE_all_strokes.txt.gz" \
    --out "${OUT_DIR}/stroke" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/hm3_no_MHC.list.txt" \
    --chunksize 500000

# -----------------------------------------------------------------------------
# Verify outputs
# -----------------------------------------------------------------------------
log_message "Verifying munged summary statistics..."
echo ""
echo "=== Munged Summary Statistics ==="
for disease in ad scz mdd bip pd ms ra asthma cad htn stroke; do
    if [[ -f "${OUT_DIR}/${disease}.sumstats.gz" ]]; then
        n_snps=$(zcat "${OUT_DIR}/${disease}.sumstats.gz" | wc -l)
        echo "[OK] ${disease}: $((n_snps - 1)) SNPs"
    else
        echo "[MISSING] ${disease}"
    fi
done

log_message "**** Job ends ****"
