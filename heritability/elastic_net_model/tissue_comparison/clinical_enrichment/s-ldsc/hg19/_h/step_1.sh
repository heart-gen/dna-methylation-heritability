#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --job-name=munge_sumstats
#SBATCH --output=logs/01.output_%j.log
#SBATCH --error=logs/01.error_%j.log

# =============================================================================
# Step 1: Munge GWAS Summary Statistics
# =============================================================================
# This script processes GWAS summary statistics for all diseases using
# LDSC's munge_sumstats.py to prepare them for S-LDSC analysis.
# =============================================================================

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

# Create output directory
OUT_DIR="./sumstats"
mkdir -p "$OUT_DIR"

# LDSC munge script
MUNGE_SCRIPT="${LDSC_DIR}/munge_sumstats.py"

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
python "$MUNGE_SCRIPT" \
    --sumstats "${GWAS_BASE}/PGC/AD/data/PGCALZ2sumstatsExcluding23andMe.txt.gz" \
    --out "${OUT_DIR}/ad" \
    --a1 testedAllele \
    --a2 otherAllele \
    --signed-sumstats z,0 \
    --p p \
    --N-col N \
    --merge-alleles "${RESOURCE_DIR}/w_hm3.snplist" \
    --chunksize 500000

# --- Schizophrenia (SCZ) ---
# Source: CLOZUK + PGC2 (Pardiñas et al. 2018)
# Columns: SNP, Freq.A1, CHR, BP, A1, A2, OR, SE, P
log_message "Processing: Schizophrenia (scz)"
python "$MUNGE_SCRIPT" \
    --sumstats "${GWAS_BASE}/PGC/SCZ/CLOZUK_PGC2/primary.qc1_filt" \
    --out "${OUT_DIR}/scz" \
    --a1 A1 \
    --a2 A2 \
    --signed-sumstats OR,1 \
    --p P \
    --snp SNP \
    --frq Freq.A1 \
    --merge-alleles "${RESOURCE_DIR}/w_hm3.snplist" \
    --chunksize 500000

# --- Major Depressive Disorder (MDD) ---
# Source: PGC MDD3 (Giannakopoulou et al. 2021)
# Columns: MarkerName, chr, pos, Allele1, Allele2, Freq1, Effect, StdErr, P.SE
log_message "Processing: Major Depressive Disorder (mdd)"
python "$MUNGE_SCRIPT" \
    --sumstats "${GWAS_BASE}/mdd/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe.txt.gz" \
    --out "${OUT_DIR}/mdd" \
    --a1 Allele1 \
    --a2 Allele2 \
    --signed-sumstats Effect,0 \
    --p P.SE \
    --snp MarkerName \
    --frq Freq1 \
    --merge-alleles "${RESOURCE_DIR}/w_hm3.snplist" \
    --chunksize 500000

# --- Bipolar Disorder (BIP) ---
# Source: PGC BIP 2021 (Mullins et al. 2021)
# Columns: #CHROM, POS, ID, A1, A2, BETA, SE, PVAL, NGT, FCAS, FCON, IMPINFO, NEFFDIV2, NCAS, NCON
log_message "Processing: Bipolar Disorder (bip)"
python "$MUNGE_SCRIPT" \
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
    --merge-alleles "${RESOURCE_DIR}/w_hm3.snplist" \
    --chunksize 500000

# --- Parkinson's Disease (PD) ---
# Source: UK Biobank imputed (self-reported)
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Parkinson's Disease (pd)"
python "$MUNGE_SCRIPT" \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_UKB_20002_1262_self_reported_parkinsons_disease.txt.gz" \
    --out "${OUT_DIR}/pd" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/w_hm3.snplist" \
    --chunksize 500000

# --- Multiple Sclerosis (MS) ---
# Source: IMMUNOBASE
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Multiple Sclerosis (ms)"
python "$MUNGE_SCRIPT" \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_IMMUNOBASE_Multiple_sclerosis_hg19.txt.gz" \
    --out "${OUT_DIR}/ms" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/w_hm3.snplist" \
    --chunksize 500000

# --- Rheumatoid Arthritis (RA) ---
# Source: OKADA trans-ethnic
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Rheumatoid Arthritis (ra)"
python "$MUNGE_SCRIPT" \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_RA_OKADA_TRANS_ETHNIC.txt.gz" \
    --out "${OUT_DIR}/ra" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/w_hm3.snplist" \
    --chunksize 500000

# --- Asthma ---
# Source: GABRIEL consortium
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Asthma (asthma)"
python "$MUNGE_SCRIPT" \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_GABRIEL_Asthma.txt.gz" \
    --out "${OUT_DIR}/asthma" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/w_hm3.snplist" \
    --chunksize 500000

# --- Coronary Artery Disease (CAD) ---
# Source: CARDIoGRAM C4D
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Coronary Artery Disease (cad)"
python "$MUNGE_SCRIPT" \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.txt.gz" \
    --out "${OUT_DIR}/cad" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/w_hm3.snplist" \
    --chunksize 500000

# --- Hypertension (HTN) ---
# Source: ICBP Systolic Blood Pressure
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Hypertension (htn)"
python "$MUNGE_SCRIPT" \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_ICBP_SystolicPressure.txt.gz" \
    --out "${OUT_DIR}/htn" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "${RESOURCE_DIR}/w_hm3.snplist" \
    --chunksize 500000

# -----------------------------------------------------------------------------
# Verify outputs
# -----------------------------------------------------------------------------
log_message "Verifying munged summary statistics..."
echo ""
echo "=== Munged Summary Statistics ==="
for disease in ad scz mdd bip pd ms ra asthma cad htn; do
    if [[ -f "${OUT_DIR}/${disease}.sumstats.gz" ]]; then
        n_snps=$(zcat "${OUT_DIR}/${disease}.sumstats.gz" | wc -l)
        echo "[OK] ${disease}: $((n_snps - 1)) SNPs"
    else
        echo "[MISSING] ${disease}"
    fi
done

log_message "**** Job ends ****"
