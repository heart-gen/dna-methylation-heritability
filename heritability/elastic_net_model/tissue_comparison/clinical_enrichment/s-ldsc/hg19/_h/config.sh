#!/bin/bash
# =============================================================================
# S-LDSC Pipeline Configuration for PSC Bridges-2
# =============================================================================
# Shared configuration file with all paths, disease/region arrays, and
# validation functions for the VMR clinical enrichment S-LDSC analysis.
# =============================================================================

# -----------------------------------------------------------------------------
# Base Paths
# -----------------------------------------------------------------------------
export PROJECT_BASE="/ocean/projects/bio250020p/kbenjamin/projects/dna-methylation-heritability"
export LDSC_DIR="/ocean/projects/bio250020p/shared/opt/ldsc"
export RESOURCE_DIR="/ocean/projects/bio250020p/shared/resources/ldsc"
export GWAS_DIR="/ocean/projects/bio250020p/shared/resources/gwas"
export LIFTOVER_DIR="/ocean/projects/bio250020p/shared/resources/liftover"

# -----------------------------------------------------------------------------
# LDSC Resource Paths
# -----------------------------------------------------------------------------
export BIM_DIR="${RESOURCE_DIR}/1000G_EUR_Phase3_plink"
export BASELINE_LD_DIR="${RESOURCE_DIR}/1000G_Phase3_baselineLD_v2.2_ldscores"
export WEIGHTS_DIR="${RESOURCE_DIR}/1000G_Phase3_weights_hm3_no_MHC"
export FRQ_DIR="${RESOURCE_DIR}/1000G_Phase3_frq"
export HAPMAP3_SNPS="${RESOURCE_DIR}/hm3_no_MHC.list.txt"

# Liftover chain file (existing file in project)
export CHAIN_FILE="/ocean/projects/bio250020p/kbenjamin/projects/dna-methylation-heritability/inputs/supportfiles/_m/hg38ToHg19.over.chain"
# Alternative if downloaded:
# export CHAIN_FILE="${LIFTOVER_DIR}/hg38ToHg19.over.chain.gz"

# -----------------------------------------------------------------------------
# Brain Regions and Heritability Status
# -----------------------------------------------------------------------------
export BRAIN_REGIONS=("caudate" "dlpfc" "hippocampus")
export HERITABILITY=("heritable_hg19" "non_heritable_hg19" "low_prediction_hg19")

# -----------------------------------------------------------------------------
# Disease/Trait Configuration (10 traits)
# -----------------------------------------------------------------------------
# Core diseases for S-LDSC analysis spanning neuronal, immune, and vascular categories

export DISEASES=("ad" "scz" "mdd" "bip" "pd" "ms" "ra" "asthma" "cad" "htn")

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
    ["htn"]="Hypertension"
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
    ["htn"]="vascular"
)
export DISEASE_CATEGORIES

# -----------------------------------------------------------------------------
# GWAS File Paths
# -----------------------------------------------------------------------------
# Note: Each GWAS file may have different column formats requiring specific
# munge_sumstats.py arguments. See step_1.sh for specific formatting.

declare -A GWAS_FILES=(
    ["ad"]="${GWAS_DIR}/PGC/AD/data/PGCALZ2sumstatsExcluding23andMe.txt.gz"
    ["scz"]="${GWAS_DIR}/PGC/SCZ/CLOZUK_PGC2/primary.qc1_filt"
    ["mdd"]="${GWAS_DIR}/mdd/jamapsy_Giannakopoulou_2021_exclude_whi_23andMe.txt.gz"
    ["bip"]="${GWAS_DIR}/bip/pgc-bip2021-all.vcf.tsv.gz"
    ["pd"]="${GWAS_DIR}/imputed_gwas_hg38_1.1/imputed_UKB_20002_1262_self_reported_parkinsons_disease.txt.gz"
    ["ms"]="${GWAS_DIR}/imputed_gwas_hg38_1.1/imputed_IMMUNOBASE_Multiple_sclerosis_hg19.txt.gz"
    ["ra"]="${GWAS_DIR}/imputed_gwas_hg38_1.1/imputed_RA_OKADA_TRANS_ETHNIC.txt.gz"
    ["asthma"]="${GWAS_DIR}/imputed_gwas_hg38_1.1/imputed_GABRIEL_Asthma.txt.gz"
    ["cad"]="${GWAS_DIR}/imputed_gwas_hg38_1.1/imputed_CARDIoGRAM_C4D_CAD_ADDITIVE.txt.gz"
    ["htn"]="${GWAS_DIR}/imputed_gwas_hg38_1.1/imputed_ICBP_SystolicPressure.txt.gz"
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
