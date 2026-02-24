#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --mem=20gb
#SBATCH --job-name=munge_sumstats
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --output=logs/munge_stats.%j.log

# =============================================================================
# Step 1: Munge GWAS Summary Statistics
# =============================================================================
# This script processes GWAS summary statistics for all diseases using
# LDSC's munge_sumstats.py to prepare them for S-LDSC analysis.
# =============================================================================

# Source configuration
SCRIPT_DIR="/projects/b1213/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/s-ldsc/hg19/_h"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

module purge
module list

log_message "**** Loading conda environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

# Create output directory
OUT_DIR="./sumstats"
mkdir -p "$OUT_DIR"

# LDSC wrapper (patches pandas/Python 3 compatibility issues)
LDSC_WRAPPER="${SCRIPT_DIR}/ldsc_wrapper.py"

# Local munge_sumstats.py (Python 3 / pandas 2.x compatible)
LOCAL_MUNGE="${SCRIPT_DIR}/munge_sumstats.py"

# GWAS directory base
GWAS_BASE="${GWAS_DIR}"

# HapMap3 SNP list with alleles (SNP, A1, A2 columns required for --merge-alleles)
HM3_SNPLIST="${RESOURCE_DIR}/w_hm3.snplist"

run_munge() {
    if ! "$@"; then
        log_message "ERROR: munge_sumstats failed"
        exit 1
    fi
}

# -----------------------------------------------------------------------------
# Munge each GWAS file
# -----------------------------------------------------------------------------
# Each GWAS has different column formats, so we handle them individually

# --- Alzheimer's Disease (AD) ---
# Source: Bellenguez et al. 2022 (GWAS Catalog harmonized file, GRCh38)
# Build handling: use harmonized rsIDs (hm_rsid), which are build-independent for LDSC.
log_message "Processing: Alzheimer's Disease (ad)"
AD_GWAS="${GWAS_BASE}/alz/bellenguez2022/35379992-GCST90027158-MONDO_0004975.h.tsv.gz"
if [[ -f "$AD_GWAS" ]]; then
    # The AD file includes both hm_x and x (duplicate columns).
    # Create a cleaned, single-rsid version to avoid ambiguous SNP column detection.
    AD_GWAS_CLEAN="${OUT_DIR}/ad.single_rsid.tsv.gz"
    if [[ ! -f "$AD_GWAS_CLEAN" ]]; then
        log_message "Cleaning AD GWAS: dropping x to keep hm_x only"
        AD_GWAS="$AD_GWAS" AD_GWAS_CLEAN="$AD_GWAS_CLEAN" python - <<'PY'
import gzip
import os

src = os.environ["AD_GWAS"]
dst = os.environ["AD_GWAS_CLEAN"]

# Define duplicate column pairs: (keep, drop)
duplicate_pairs = [
    ("hm_rsid", "variant_id"),
    ("hm_chrom", "chromosome"),
    ("hm_pos", "base_pair_location"),
    ("hm_other_allele", "other_allele"),
    ("hm_effect_allele", "effect_allele"),
    ("hm_beta", "beta"),
    ("hm_odds_ratio", "odds_ratio"),
    ("hm_ci_lower", "ci_lower"),
    ("hm_ci_upper", "ci_upper"),
    ("hm_effect_allele_frequency", "effect_allele_frequency"),
]

with gzip.open(src, "rt") as fin, gzip.open(dst, "wt") as fout:
    header = fin.readline().rstrip("\n").split("\t")

    # Determine which columns to drop
    drop_cols = set()
    for keep, drop in duplicate_pairs:
        if keep in header and drop in header:
            drop_cols.add(drop)

    keep_idx = [i for i, col in enumerate(header) if col not in drop_cols]

    # Write cleaned header
    fout.write("\t".join(header[i] for i in keep_idx) + "\n")

    # Write cleaned rows
    for line in fin:
        parts = line.rstrip("\n").split("\t")
        fout.write("\t".join(parts[i] for i in keep_idx) + "\n")
PY
    fi
    run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
        --sumstats "$AD_GWAS_CLEAN" \
        --out "${OUT_DIR}/ad" \
        --a1 hm_effect_allele \
        --a2 hm_other_allele \
        --signed-sumstats hm_beta,0 \
        --p p_value \
        --snp hm_rsid \
        --frq hm_effect_allele_frequency \
        --N-cas-col n_cas \
        --N-con-col n_con \
        --merge-alleles "$HM3_SNPLIST" \
        --chunksize 500000
else
    log_message "WARNING: AD GWAS not found (${AD_GWAS}). Skipping until download completes."
fi

# --- Schizophrenia (SCZ) ---
# Source: PGC3 (Trubetskoy et al. 2022 Nature)
# Columns: CHROM, ID, POS, A1, A2, FCAS, FCON, IMPINFO, BETA, SE, PVAL, NCAS, NCON, NEFF
# NOTE: VCF format with 73 ## header lines that must be skipped
log_message "Processing: Schizophrenia (scz)"
run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
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
run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
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
run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
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
# Source: GCST009325 harmonized file (GRCh38)
# Build handling: use rsid column to avoid coordinate-build mismatch in hg19 LDSC pipeline.
log_message "Processing: Parkinson's Disease (pd)"
PD_GWAS="${GWAS_BASE}/PD/data/GCST009325.h.tsv.gz"
if [[ -f "$PD_GWAS" ]]; then
    # The PD file includes both hm_x and x (duplicate columns).
    # Create a cleaned, single-rsid version to avoid ambiguous SNP column detection.
    PD_GWAS_CLEAN="${OUT_DIR}/pd.single_rsid.tsv.gz"
    if [[ ! -f "$PD_GWAS_CLEAN" ]]; then
        log_message "Cleaning PD GWAS: dropping x to keep hm_x only"
        PD_GWAS="$PD_GWAS" PD_GWAS_CLEAN="$PD_GWAS_CLEAN" python - <<'PY'
import gzip
import os

src = os.environ["PD_GWAS"]
dst = os.environ["PD_GWAS_CLEAN"]

# Define duplicate column pairs: (keep, drop)
duplicate_pairs = [
    ("rsid", "variant_id"),
]

with gzip.open(src, "rt") as fin, gzip.open(dst, "wt") as fout:
    header = fin.readline().rstrip("\n").split("\t")

    # Determine which columns to drop
    drop_cols = set()
    for keep, drop in duplicate_pairs:
        if keep in header and drop in header:
            drop_cols.add(drop)

    keep_idx = [i for i, col in enumerate(header) if col not in drop_cols]

    # Write cleaned header
    fout.write("\t".join(header[i] for i in keep_idx) + "\n")

    # Write cleaned rows
    for line in fin:
        parts = line.rstrip("\n").split("\t")
        fout.write("\t".join(parts[i] for i in keep_idx) + "\n")
PY
    fi
    run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
        --sumstats "$PD_GWAS_CLEAN" \
        --out "${OUT_DIR}/pd" \
        --a1 effect_allele \
        --a2 other_allele \
        --signed-sumstats beta,0 \
        --p p_value \
        --snp rsid \
        --frq effect_allele_frequency \
        --N-cas-col N_cases \
        --N-con-col N_controls \
        --merge-alleles "$HM3_SNPLIST" \
        --chunksize 500000
else
    log_message "WARNING: PD GWAS not found (${PD_GWAS}). Skipping until download completes."
fi

# --- Multiple Sclerosis (MS) ---
# Source: IMMUNOBASE
# Columns: variant_id, panel_variant_id, chromosome, position, effect_allele, non_effect_allele, current_build, frequency, sample_size, zscore, pvalue
log_message "Processing: Multiple Sclerosis (ms)"
run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
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
run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
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
run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
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
run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
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
# Source: Malik et al. 2018 ischemic stroke (build37 formatted GWAS Catalog file)
# Build handling: this specific file is build37/hg19, so no liftover is required.
log_message "Processing: Stroke (stroke)"
STROKE_GWAS="${GWAS_BASE}/stroke/29531354-GCST005843-HP_0002140-build37.f.tsv.gz"
if [[ -f "$STROKE_GWAS" ]]; then
    run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
        --sumstats "$STROKE_GWAS" \
        --out "${OUT_DIR}/stroke" \
        --a1 effect_allele \
        --a2 other_allele \
        --signed-sumstats beta,0 \
        --p p_value \
        --snp variant_id \
        --frq effect_allele_frequency \
        --N 517525 \
        --merge-alleles "$HM3_SNPLIST" \
        --chunksize 500000
else
    log_message "WARNING: Stroke GWAS not found (${STROKE_GWAS}). Skipping until download completes."
fi

# --- Height (Control) ---
# Source: UK Biobank standing height (hg38 imputed file)
# No HM3 merge.
log_message "Processing: Height (height control)"
HEIGHT_GWAS="${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_UKB_50_Standing_height.txt.gz"
CHAIN_FILE=/projects/b1213/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/liftOver/hg38ToHg19.over.chain.gz

zcat "$HEIGHT_GWAS" | \
awk 'NR>1 {print $3"\t"$4-1"\t"$4"\t"$1}' \
> $OUT_DIR/height.bed

liftOver $OUT_DIR/height.bed $CHAIN_FILE $OUT_DIR/height_lifted.bed $OUT_DIR/height_unmapped.bed

awk '{print $4"\t"$3}' $OUT_DIR/height_lifted.bed > $OUT_DIR/height_new_positions.txt

zcat "$HEIGHT_GWAS" | \
awk 'BEGIN{OFS="\t"}
NR==FNR {pos[$1]=$2; next}
NR==1 {print; next}
{
    if($1 in pos) $4=pos[$1];
    print
}' $OUT_DIR/height_new_positions.txt - | \
gzip > $OUT_DIR/height_lifted.txt.gz

if [[ -f "$OUT_DIR/height_lifted.txt.gz" ]]; then
    run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
        --sumstats "$OUT_DIR/height_lifted.txt.gz" \
        --out "${OUT_DIR}/height" \
        --a1 effect_allele \
        --a2 non_effect_allele \
        --signed-sumstats zscore,0 \
        --p pvalue \
        --snp variant_id \
        --frq frequency \
        --N-col sample_size \
        --merge-alleles "$HM3_SNPLIST" \
        --chunksize 500000 \
else
    log_message "WARNING: Height GWAS not found (${OUT_DIR}/height_lifted.txt.gz). Skipping."
fi

# -----------------------------------------------------------------------------
# Verify outputs
# -----------------------------------------------------------------------------
log_message "Verifying munged summary statistics..."
echo ""
echo "=== Munged Summary Statistics ==="
for disease in ad scz mdd bip pd ms ra asthma cad stroke height; do
    if [[ -f "${OUT_DIR}/${disease}.sumstats.gz" ]]; then
        n_snps=$(zcat "${OUT_DIR}/${disease}.sumstats.gz" | wc -l)
        echo "[OK] ${disease}: $((n_snps - 1)) SNPs"
    else
        echo "[MISSING] ${disease}"
    fi
done

conda deactivate
log_message "**** Job ends ****"