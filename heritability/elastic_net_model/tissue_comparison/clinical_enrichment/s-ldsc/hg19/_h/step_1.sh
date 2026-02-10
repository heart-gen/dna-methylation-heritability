#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=20gb
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
PROJECT_BASE="/path/to/dna-methylation-heritability"
SCRIPT_DIR="${PROJECT_BASE}/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/s-ldsc/hg19/_h"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading conda environment ****"
conda activate /projects/p32505/opt/envs/genomics

# Create output directory
OUT_DIR="./sumstats"
mkdir -p "$OUT_DIR"

# LDSC wrapper (patches pandas/Python 3 compatibility issues)
LDSC_WRAPPER="${SCRIPT_DIR}/ldsc_wrapper.py"

# Local munge_sumstats.py (Python 3 / pandas 2.x compatible)
LOCAL_MUNGE="${SCRIPT_DIR}/munge_sumstats.py"

# GWAS directory base
GWAS_BASE="projects/b1213/resources/gwas"

# HapMap3 SNP list with alleles (SNP, A1, A2 columns required for --merge-alleles)
HM3_SNPLIST="${RESOURCE_DIR}/w_hm3.snplist"

# -----------------------------------------------------------------------------
# Munge each GWAS file
# -----------------------------------------------------------------------------
# Each GWAS has different column formats, so we handle them individually

# --- Alzheimer's Disease (AD) ---
# Source: PGC ALZ2 (Kunkle et al. 2019)
# Columns: chr, PosGRCh37, testedAllele, otherAllele, z, p, N
# NOTE: This file has NO SNP column (rs IDs) - it only has chr:position
# We need to use a pre-processed file or skip this trait
# For now, using a different Alzheimer's GWAS that has rs IDs
log_message "Processing: Alzheimer's Disease (ad)"
# Check if we have a sumstats file with SNP IDs
if [[ -f "${GWAS_BASE}/PGC/AD/data/PGCALZ2sumstatsExcluding23andMe.with_rsid.txt.gz" ]]; then
    python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
        --sumstats "${GWAS_BASE}/PGC/AD/data/PGCALZ2sumstatsExcluding23andMe.with_rsid.txt.gz" \
        --out "${OUT_DIR}/ad" \
        --a1 testedAllele \
        --a2 otherAllele \
        --signed-sumstats z,0 \
        --p p \
        --N-col N \
        --merge-alleles "$HM3_SNPLIST" \
        --chunksize 500000
else
    # Use the pre-munged PASS Alzheimer's sumstats if available
    if [[ -f "${RESOURCE_DIR}/sumstats/PASS_Alzheimer.sumstats.gz" ]]; then
        log_message "Using pre-munged PASS_Alzheimer.sumstats.gz"
        cp "${RESOURCE_DIR}/sumstats/PASS_Alzheimer.sumstats.gz" "${OUT_DIR}/ad.sumstats.gz"
    else
        log_message "WARNING: AD GWAS file lacks SNP column and no alternative found. Skipping."
    fi
fi

# --- Schizophrenia (SCZ) ---
# Source: PGC3 (Trubetskoy et al. 2022 Nature)
# Columns: CHROM, ID, POS, A1, A2, FCAS, FCON, IMPINFO, BETA, SE, PVAL, NCAS, NCON, NEFF
# NOTE: VCF format with 73 ## header lines that must be skipped
log_message "Processing: Schizophrenia (scz)"
python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
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
    --merge-alleles "$HM3_SNPLIST" \
    --chunksize 500000 \
    --skip-rows 73

# --- Major Depressive Disorder (MDD) ---
# Source: PGC MDD3 (Giannakopoulou et al. 2021)
# Columns: MarkerName, chr, pos, Allele1, Allele2, Freq1, Effect, StdErr, P.SE
# NOTE: File lacks N column. Using N from paper: 166,773 cases + 507,679 controls
log_message "Processing: Major Depressive Disorder (mdd)"
python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
    --sumstats "${GWAS_BASE}/mdd/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe.txt.gz" \
    --out "${OUT_DIR}/mdd" \
    --a1 Allele1 \
    --a2 Allele2 \
    --signed-sumstats Effect,0 \
    --p P.SE \
    --snp MarkerName \
    --frq Freq1 \
    --N-cas 166773 \
    --N-con 507679 \
    --merge-alleles "$HM3_SNPLIST" \
    --chunksize 500000

# --- Bipolar Disorder (BIP) ---
# Source: PGC BIP 2021 (Mullins et al. 2021)
# Columns: #CHROM, POS, ID, A1, A2, BETA, SE, PVAL, NGT, FCAS, FCON, IMPINFO, NEFFDIV2, NCAS, NCON
# NOTE: VCF format with 72 ## header lines that must be skipped
log_message "Processing: Bipolar Disorder (bip)"
python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
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
    --merge-alleles "$HM3_SNPLIST" \
    --chunksize 500000 \
    --skip-rows 72

# --- Parkinson's Disease (PD) ---
# Source: UK Biobank imputed (self-reported)
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
# NOTE: File is in hg38 coordinates - requires liftover to hg19 before munging
log_message "Processing: Parkinson's Disease (pd)"
PD_GWAS_HG38="${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_UKB_20002_1262_self_reported_parkinsons_disease.txt.gz"
PD_GWAS_HG19="${OUT_DIR}/pd_hg19.txt.gz"
CHAIN_FILE="${SCRIPT_DIR}/../../../../../../../inputs/supportfiles/_m/hg38ToHg19.over.chain"

# Liftover PD GWAS from hg38 to hg19
if [[ ! -f "$PD_GWAS_HG19" ]]; then
    log_message "  Lifting over PD GWAS from hg38 to hg19..."
    python -c "
import pandas as pd
import gzip
from pyliftover import LiftOver

lo = LiftOver('${CHAIN_FILE}')

def liftover_pos(chrom, pos):
    try:
        chrom_str = f'chr{chrom}' if not str(chrom).startswith('chr') else str(chrom)
        result = lo.convert_coordinate(chrom_str, int(pos))
        if result and len(result) > 0:
            return int(result[0][1])
        return None
    except:
        return None

# Read in chunks to handle large file
chunks = []
for chunk in pd.read_csv('${PD_GWAS_HG38}', sep='\t', chunksize=500000):
    chunk['position_hg19'] = chunk.apply(lambda x: liftover_pos(x['chromosome'], x['position']), axis=1)
    chunk = chunk.dropna(subset=['position_hg19'])
    chunk['position'] = chunk['position_hg19'].astype(int)
    chunk = chunk.drop(columns=['position_hg19'])
    chunks.append(chunk)
    print(f'  Processed chunk: {len(chunk)} SNPs retained')

df = pd.concat(chunks, ignore_index=True)
print(f'  Total SNPs after liftover: {len(df)}')
df.to_csv('${PD_GWAS_HG19}', sep='\t', index=False, compression='gzip', na_rep='NA')
print('  Liftover complete!')
"
fi

python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
    --sumstats "$PD_GWAS_HG19" \
    --out "${OUT_DIR}/pd" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "$HM3_SNPLIST" \
    --chunksize 500000

# --- Multiple Sclerosis (MS) ---
# Source: IMMUNOBASE
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Multiple Sclerosis (ms)"
python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_IMMUNOBASE_Multiple_sclerosis_hg19.txt.gz" \
    --out "${OUT_DIR}/ms" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "$HM3_SNPLIST" \
    --chunksize 500000

# --- Rheumatoid Arthritis (RA) ---
# Source: OKADA trans-ethnic
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Rheumatoid Arthritis (ra)"
python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_RA_OKADA_TRANS_ETHNIC.txt.gz" \
    --out "${OUT_DIR}/ra" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "$HM3_SNPLIST" \
    --chunksize 500000

# --- Asthma ---
# Source: GABRIEL consortium
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Asthma (asthma)"
python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_GABRIEL_Asthma.txt.gz" \
    --out "${OUT_DIR}/asthma" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "$HM3_SNPLIST" \
    --chunksize 500000

# --- Coronary Artery Disease (CAD) ---
# Source: CARDIoGRAM C4D
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Coronary Artery Disease (cad)"
python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.txt.gz" \
    --out "${OUT_DIR}/cad" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "$HM3_SNPLIST" \
    --chunksize 500000

# --- Hypertension (HTN) ---
# REMOVED: No hg19 GWAS available. The imputed_ICBP_SystolicPressure.txt.gz is hg38
# and has no matching SNPs in the hg19 HapMap3 reference.
# Replaced with Stroke as the second vascular trait.

# --- Stroke ---
# Source: ISGC METASTROKE (Malik et al. 2016) - All strokes
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Stroke (stroke)"
python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
    --sumstats "${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_ISGC_Malik_2016_METASTROKE_all_strokes.txt.gz" \
    --out "${OUT_DIR}/stroke" \
    --a1 effect_allele \
    --a2 non_effect_allele \
    --signed-sumstats zscore,0 \
    --p pvalue \
    --snp variant_id \
    --frq frequency \
    --N-col sample_size \
    --merge-alleles "$HM3_SNPLIST" \
    --chunksize 500000

# -----------------------------------------------------------------------------
# Verify outputs
# -----------------------------------------------------------------------------
log_message "Verifying munged summary statistics..."
echo ""
echo "=== Munged Summary Statistics ==="
# Note: htn replaced with stroke as vascular trait (no hg19 HTN GWAS available)
for disease in ad scz mdd bip pd ms ra asthma cad stroke; do
    if [[ -f "${OUT_DIR}/${disease}.sumstats.gz" ]]; then
        n_snps=$(zcat "${OUT_DIR}/${disease}.sumstats.gz" | wc -l)
        echo "[OK] ${disease}: $((n_snps - 1)) SNPs"
    else
        echo "[MISSING] ${disease}"
    fi
done

log_message "**** Job ends ****"
