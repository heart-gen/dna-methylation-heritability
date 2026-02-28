#!/bin/bash
# =============================================================================
# S-LDSC Pipeline Configuration for QUEST
# =============================================================================
# Shared configuration file with all paths, disease/region arrays, and
# validation functions for the VMR clinical enrichment S-LDSC analysis.
# =============================================================================

# -----------------------------------------------------------------------------
# Base Paths
# -----------------------------------------------------------------------------
CONFIG_DIR="/projects/b1213/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/clinical_enrichment/s-ldsc/hg19/_h"
export PROJECT_BASE="${PROJECT_BASE:-$(cd "${CONFIG_DIR}/../../../../../../.." && pwd)}"
export LDSC_DIR="/projects/p32505/opt/ldsc"
export RESOURCE_DIR="/projects/b1213/resources/ldsc"
export GWAS_DIR="/projects/b1213/resources/gwas"

# -----------------------------------------------------------------------------
# LDSC Resource Paths
# -----------------------------------------------------------------------------
export BIM_DIR="${RESOURCE_DIR}/1000G_EUR_Phase3_plink"
export BASELINE_LD_DIR="${RESOURCE_DIR}/1000G_Phase3_baselineLD_v2.2_ldscores"
export WEIGHTS_DIR="${RESOURCE_DIR}/1000G_Phase3_weights_hm3_no_MHC"
export FRQ_DIR="${RESOURCE_DIR}/1000G_Phase3_frq"
# HapMap3 SNP list with alleles (required for --merge-alleles in munge_sumstats)
# The w_hm3.snplist file has SNP, A1, A2 columns
export HAPMAP3_SNPS_ALLELES="${RESOURCE_DIR}/w_hm3.snplist"
# HapMap3 SNP list with just SNP IDs (required for --print-snps in ldsc.py)
# Single-column format, no header
export HAPMAP3_SNPS="${RESOURCE_DIR}/hm3_no_MHC.list.txt"

# Liftover chain file (existing file in project)
export CHAIN_FILE="${PROJECT_BASE}/inputs/supportfiles/_m/hg38ToHg19.over.chain"

# -----------------------------------------------------------------------------
# Brain Regions and Heritability Status
# -----------------------------------------------------------------------------
export BRAIN_REGIONS=("caudate" "dlpfc" "hippocampus")
export HERITABILITY=("heritable_hg19" "non_heritable_hg19" "low_prediction_hg19")

# -----------------------------------------------------------------------------
# Disease/Trait Configuration (11 traits)
# -----------------------------------------------------------------------------
# Core diseases for S-LDSC analysis spanning neuronal, immune, vascular,
# and control categories.

export DISEASES=("ad" "scz" "mdd" "bip" "pd" "ms" "ra" "asthma" "cad" "stroke" "smoking" "height")

# Full disease names for reporting
declare -A DISEASE_NAMES=(
    ["ad"]="Alzheimer's Disease"
    ["scz"]="Schizophrenia"
    ["mdd"]="Major Depressive Disorder"
    ["bip"]="Bipolar Disorder"
    ["pd"]="Parkinson's Disease"
    ["ms"]="Multiple Sclerosis"
    ["ra"]="Rheumatoid Arthritis"
    ["asthma"]="Asthma"
    ["cad"]="Coronary Artery Disease"
    ["stroke"]="Stroke"
    ["substance_abuse"]="Substance Abuse"
    ["smoking"]="Ever Smoker"
    ["height"]="Standing Height (Control)"
)
export DISEASE_NAMES

# Disease categories for grouping
declare -A DISEASE_CATEGORIES=(
    ["ad"]="neuronal"
    ["scz"]="neuronal"
    ["mdd"]="neuronal"
    ["bip"]="neuronal"
    ["pd"]="neuronal"
    ["ms"]="immune"
    ["ra"]="immune"
    ["asthma"]="immune"
    ["cad"]="vascular"
    ["stroke"]="vascular"
    ["smoking"]="tbd"
    ["substance_abuse"]="tbd"
    ["height"]="control"
)
export DISEASE_CATEGORIES

# -----------------------------------------------------------------------------
# GWAS File Paths
# -----------------------------------------------------------------------------
# Note: Each GWAS file may have different column formats requiring specific
# munge_sumstats.py arguments. See step_1.sh for specific formatting.

declare -A GWAS_FILES=(
    ["ad"]="${GWAS_DIR}/alz/bellenguez2022/35379992-GCST90027158-MONDO_0004975.h.tsv.gz"
    ["scz"]="${GWAS_DIR}/PGC/SCZ/PGC3/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
    ["mdd"]="${GWAS_DIR}/mdd/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe.txt.gz"
    ["bip"]="${GWAS_DIR}/bip/pgc-bip2021-all.vcf.tsv.gz"
    ["pd"]="${GWAS_DIR}/PD/data/GCST009325.h.tsv.gz"
    ["ms"]="${GWAS_DIR}/imputed_gwas_hg38_1.1/imputed_IMMUNOBASE_Multiple_sclerosis_hg19.txt.gz"
    ["ra"]="${GWAS_DIR}/imputed_gwas_hg38_1.1/imputed_RA_OKADA_TRANS_ETHNIC.txt.gz"
    ["asthma"]="${GWAS_DIR}/imputed_gwas_hg38_1.1/imputed_GABRIEL_Asthma.txt.gz"
    ["cad"]="${GWAS_DIR}/imputed_gwas_hg38_1.1/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.txt.gz"
    ["stroke"]="${GWAS_DIR}/stroke/29531354-GCST005843-HP_0002140-build37.f.tsv.gz"
    ["smoking"]="../EVER_SMOKER_GWAS_MA_UKB+TAG.txt"
    ["substance_abuse"]="../hg19/GCST90435891.tsv.gz"
    ["height"]="${GWAS_DIR}/imputed_gwas_hg38_1.1/imputed_UKB_50_Standing_height.txt.gz"
)
export GWAS_FILES

# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

# Log function with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}
export -f log_message

# Print SLURM job info
print_job_info() {
    echo "**** PSC Bridges-2 info ****"
    echo "User: ${USER}"
    echo "Job id: ${SLURM_JOBID}"
    echo "Job name: ${SLURM_JOB_NAME}"
    echo "Node name: ${SLURM_NODENAME}"
    echo "Hostname: ${HOSTNAME}"
    echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"
}
export -f print_job_info

# Validate that required resource directories exist
validate_resources() {
    local missing=0

    echo "Validating resources..."

    if [[ ! -d "$LDSC_DIR" ]]; then
        echo "[ERROR] LDSC directory not found: $LDSC_DIR"
        missing=1
    fi

    if [[ ! -d "$BIM_DIR" ]]; then
        echo "[ERROR] BIM directory not found: $BIM_DIR"
        echo "       Run download_resources.sh first"
        missing=1
    fi

    if [[ ! -d "$BASELINE_LD_DIR" ]]; then
        echo "[ERROR] Baseline LD directory not found: $BASELINE_LD_DIR"
        echo "       Run download_resources.sh first"
        missing=1
    fi

    if [[ ! -d "$WEIGHTS_DIR" ]]; then
        echo "[ERROR] Weights directory not found: $WEIGHTS_DIR"
        echo "       Run download_resources.sh first"
        missing=1
    fi

    if [[ ! -d "$FRQ_DIR" ]]; then
        echo "[ERROR] Frequency directory not found: $FRQ_DIR"
        echo "       Run download_resources.sh first"
        missing=1
    fi

    if [[ ! -f "$CHAIN_FILE" ]]; then
        echo "[ERROR] Chain file not found: $CHAIN_FILE"
        missing=1
    fi

    if [[ $missing -eq 1 ]]; then
        echo ""
        echo "Some required resources are missing. Please ensure all resources are downloaded."
        return 1
    fi

    echo "All resources validated successfully."
    return 0
}
export -f validate_resources

# Get input file path for a brain region
get_input_file() {
    local region="$1"
    echo "${PROJECT_BASE}/heritability/elastic_net_model/${region}/_m/${region}_summary_elastic-net.tsv"
}
export -f get_input_file
