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

SCRIPT_DIR="/projects/b1213/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/BA_only/tissue_comparison/sensitivity_analysis/clinical_enrichment/s-ldsc/hg19/_h"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

module purge
module list

log_message "**** Loading conda environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

OUT_DIR="./sumstats"
mkdir -p "$OUT_DIR"

LDSC_WRAPPER="${SCRIPT_DIR}/ldsc_wrapper.py"
LOCAL_MUNGE="${SCRIPT_DIR}/munge_sumstats.py"
GWAS_BASE="${GWAS_DIR}"
HM3_SNPLIST="${RESOURCE_DIR}/w_hm3.snplist"

run_munge() {
    if ! "$@"; then
        log_message "ERROR: munge_sumstats failed"
        exit 1
    fi
}

# -----------------------------------------------------------------------------
# --- Alzheimer's Disease (AD) ---
log_message "Processing: Alzheimer's Disease (ad)"
AD_GWAS="${GWAS_BASE}/alz/bellenguez2022/35379992-GCST90027158-MONDO_0004975.h.tsv.gz"
if [[ -f "$AD_GWAS" ]]; then
    AD_GWAS_CLEAN="${OUT_DIR}/ad.single_rsid.tsv.gz"
    if [[ ! -f "$AD_GWAS_CLEAN" ]]; then
        log_message "Cleaning AD GWAS: dropping duplicate columns"
        AD_GWAS="$AD_GWAS" AD_GWAS_CLEAN="$AD_GWAS_CLEAN" python - <<'PY'
import gzip, os
src = os.environ["AD_GWAS"]
dst = os.environ["AD_GWAS_CLEAN"]
duplicate_pairs = [
    ("hm_rsid","variant_id"), ("hm_chrom","chromosome"), ("hm_pos","base_pair_location"),
    ("hm_other_allele","other_allele"), ("hm_effect_allele","effect_allele"),
    ("hm_beta","beta"), ("hm_odds_ratio","odds_ratio"), ("hm_ci_lower","ci_lower"),
    ("hm_ci_upper","ci_upper"), ("hm_effect_allele_frequency","effect_allele_frequency")
]
with gzip.open(src,"rt") as fin, gzip.open(dst,"wt") as fout:
    header = fin.readline().rstrip("\n").split("\t")
    drop_cols = {drop for keep,drop in duplicate_pairs if keep in header and drop in header}
    keep_idx = [i for i,col in enumerate(header) if col not in drop_cols]
    fout.write("\t".join(header[i] for i in keep_idx)+"\n")
    for line in fin:
        parts = line.rstrip("\n").split("\t")
        fout.write("\t".join(parts[i] for i in keep_idx)+"\n")
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
    log_message "WARNING: AD GWAS not found (${AD_GWAS}). Skipping."
fi

# -----------------------------------------------------------------------------
# --- Schizophrenia (SCZ) ---
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

# -----------------------------------------------------------------------------
# --- Major Depressive Disorder (MDD) ---
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

# -----------------------------------------------------------------------------
# --- Bipolar Disorder (BIP) ---
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

# -----------------------------------------------------------------------------
# --- Parkinson's Disease (PD) ---
log_message "Processing: Parkinson's Disease (pd)"
PD_GWAS="${GWAS_BASE}/PD/data/GCST009325.h.tsv.gz"
if [[ -f "$PD_GWAS" ]]; then
    PD_GWAS_CLEAN="${OUT_DIR}/pd.single_rsid.tsv.gz"
    if [[ ! -f "$PD_GWAS_CLEAN" ]]; then
        log_message "Cleaning PD GWAS: dropping variant_id"
        PD_GWAS="$PD_GWAS" PD_GWAS_CLEAN="$PD_GWAS_CLEAN" python - <<'PY'
import gzip, os
src = os.environ["PD_GWAS"]
dst = os.environ["PD_GWAS_CLEAN"]
duplicate_pairs = [("rsid","variant_id")]
with gzip.open(src,"rt") as fin, gzip.open(dst,"wt") as fout:
    header = fin.readline().rstrip("\n").split("\t")
    drop_cols = {drop for keep,drop in duplicate_pairs if keep in header and drop in header}
    keep_idx = [i for i,col in enumerate(header) if col not in drop_cols]
    fout.write("\t".join(header[i] for i in keep_idx)+"\n")
    for line in fin:
        parts = line.rstrip("\n").split("\t")
        fout.write("\t".join(parts[i] for i in keep_idx)+"\n")
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
    log_message "WARNING: PD GWAS not found (${PD_GWAS}). Skipping."
fi

# -----------------------------------------------------------------------------
# --- Height (Control) ---
log_message "Processing: Height (control)"
HEIGHT_GWAS="../GCST90468178.tsv.gz"
if [[ -f "$HEIGHT_GWAS" ]]; then
    HEIGHT_GWAS_CLEAN="${OUT_DIR}/height.single_rsid.tsv.gz"
    if [[ ! -f "$HEIGHT_GWAS_CLEAN" ]]; then
        log_message "Cleaning Height GWAS: dropping variant_id"
        HEIGHT_GWAS="$HEIGHT_GWAS" HEIGHT_GWAS_CLEAN="$HEIGHT_GWAS_CLEAN" python - <<'PY'
import gzip, os

src = os.environ["HEIGHT_GWAS"]
dst = os.environ["HEIGHT_GWAS_CLEAN"]

with gzip.open(src,"rt") as fin, gzip.open(dst,"wt") as fout:
    header = fin.readline().rstrip("\n").split("\t")

    # Drop variant_id unconditionally if present
    keep_idx = [i for i,col in enumerate(header) if col != "variant_id"]

    fout.write("\t".join(header[i] for i in keep_idx)+"\n")

    for line in fin:
        parts = line.rstrip("\n").split("\t")
        fout.write("\t".join(parts[i] for i in keep_idx)+"\n")
PY
    fi
    run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
        --sumstats "$HEIGHT_GWAS_CLEAN" \
        --out "${OUT_DIR}/height" \
        --a1 effect_allele \
        --a2 other_allele \
        --signed-sumstats beta,0 \
        --p p_value \
        --snp rs_id \
        --frq effect_allele_frequency \
        --N-col n \
        --merge-alleles "$HM3_SNPLIST" \
        --chunksize 500000
else
    log_message "WARNING: Height GWAS not found ($HEIGHT_GWAS_CLEAN). Skipping."
fi

# -----------------------------------------------------------------------------
# --- Other GWAS (MS, RA, Asthma, CAD, Stroke, Smoking, Substance Abuse) ---
# Pattern: all direct run_munge calls (no heredoc)
declare -A gwas_files=(
    ["ms"]="${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_IMMUNOBASE_Multiple_sclerosis_hg19.txt.gz"
    ["ra"]="${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_RA_OKADA_TRANS_ETHNIC.txt.gz"
    ["asthma"]="${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_GABRIEL_Asthma.txt.gz"
    ["cad"]="${GWAS_BASE}/imputed_gwas_hg38_1.1/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.txt.gz"
    ["stroke"]="${GWAS_BASE}/stroke/29531354-GCST005843-HP_0002140-build37.f.tsv.gz"
    ["smoking"]="../EVER_SMOKER_GWAS_MA_UKB+TAG.txt.gz"
    ["substance_abuse"]="../GCST90435891.tsv.gz"
)

for disease in "${!gwas_files[@]}"; do
    FILE="${gwas_files[$disease]}"
    if [[ -f "$FILE" ]]; then
        log_message "Processing: $disease"
        case "$disease" in
            "ms"|"ra"|"asthma"|"cad")
                run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
                    --sumstats "$FILE" \
                    --out "${OUT_DIR}/$disease" \
                    --a1 effect_allele \
                    --a2 non_effect_allele \
                    --signed-sumstats zscore,0 \
                    --p pvalue \
                    --snp variant_id \
                    --frq frequency \
                    --N-col sample_size \
                    --merge-alleles "$HM3_SNPLIST" \
                    --chunksize 500000
                ;;
            "stroke")
                run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
                    --sumstats "$FILE" \
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
                ;;
            "smoking")
                run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
                    --sumstats "$FILE" \
                    --out "${OUT_DIR}/smoking" \
                    --a1 A1 \
                    --a2 A2 \
                    --signed-sumstats Beta,0 \
                    --p Pval \
                    --snp MarkerName \
                    --frq EAF_A1 \
                    --N 518633 \
                    --merge-alleles "$HM3_SNPLIST" \
                    --chunksize 500000
                ;;
            "substance_abuse")
                run_munge python "$LDSC_WRAPPER" "$LDSC_DIR" "$LOCAL_MUNGE" \
                    --sumstats "$FILE" \
                    --out "${OUT_DIR}/substance_abuse" \
                    --a1 effect_allele \
                    --a2 other_allele \
                    --signed-sumstats beta,0 \
                    --p p_value \
                    --snp variant_id \
                    --frq effect_allele_frequency \
                    --N 518633 \
                    --merge-alleles "$HM3_SNPLIST" \
                    --chunksize 500000
                ;;
        esac
    else
        log_message "WARNING: $disease GWAS not found ($FILE). Skipping."
    fi
done

# -----------------------------------------------------------------------------
# Verify outputs
log_message "Verifying munged summary statistics..."
echo ""
for disease in ad scz mdd bip pd ms ra asthma cad stroke smoking substance_abuse height; do
    if [[ -f "${OUT_DIR}/${disease}.sumstats.gz" ]]; then
        n_snps=$(zcat "${OUT_DIR}/${disease}.sumstats.gz" | wc -l)
        echo "[OK] $disease: $((n_snps - 1)) SNPs"
    else
        echo "[MISSING] $disease"
    fi
done

conda deactivate
log_message "**** Job ends ****"