#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=08:00:00
#SBATCH --mem=20gb
#SBATCH --job-name=partition_heritability
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=elisajohnson2027@u.northwestern.edu
#SBATCH --output=logs/partition_heritability.%j.log

# =============================================================================
# Step 4: Partition Heritability with S-LDSC
# =============================================================================
# This script runs stratified LD score regression (S-LDSC) to partition
# heritability for each disease/trait across brain regions
# Processes 11 traits x 3 regions = 33 analyses.
# =============================================================================

# Source configuration
SCRIPT_DIR="/projects/b1213/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/all_individuals/tissue_comparison/clinical_enrichment/s-ldsc/hg19/_h"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

module purge
module list

log_message "**** Loading conda environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

# Validate resources
if ! validate_resources; then
    log_message "ERROR: Resource validation failed. Exiting."
    exit 1
fi

# LDSC wrapper (patches pandas/Python 3 compatibility issues)
LDSC_WRAPPER="${SCRIPT_DIR}/ldsc_wrapper.py"

# Resource paths
BASELINE_LD_DIR="${RESOURCE_DIR}/1000G_Phase3_baselineLD_v2.2_ldscores"
WEIGHTS_DIR="${RESOURCE_DIR}/1000G_Phase3_weights_hm3_no_MHC"
FRQ_DIR="${RESOURCE_DIR}/1000G_Phase3_frq"

# Summary statistics directory
SUMSTATS_DIR="./sumstats"

# DISEASES array is defined in config.sh

# Counter for progress
total=$((${#DISEASES[@]} * ${#BRAIN_REGIONS[@]}))
current=0

# Process each disease
for DISEASE in "${DISEASES[@]}"; do
    log_message "Processing disease: ${DISEASE}"

    # Check if sumstats file exists
    SUMSTATS_FILE="${SUMSTATS_DIR}/${DISEASE}.sumstats.gz"
    if [[ ! -f "$SUMSTATS_FILE" ]]; then
        log_message "  WARNING: Summary statistics not found: $SUMSTATS_FILE"
        log_message "           Run step_1.sh first to munge GWAS data"
        continue
    fi

    # Process each brain region
    for REGION in "${BRAIN_REGIONS[@]}"; do
    ((current++))
        log_message "  Processing region: ${REGION} (${current}/${total})"

        # Custom LD scores directory
        CUSTOM_LD_DIR="./custom_ldscores/${REGION}/EA"

        # Output directory
        OUT_DIR="./results/${DISEASE}/${REGION}/EA"
        mkdir -p "$OUT_DIR"

        # Output file prefix
        OUT_PREFIX="${OUT_DIR}/${DISEASE}_${REGION}_EA"

        # Check if custom LD scores exist
        if [[ ! -f "${CUSTOM_LD_DIR}/${REGION}.1.l2.ldscore.gz" ]]; then
            log_message "      WARNING: Custom LD scores not found in: $CUSTOM_LD_DIR"
            log_message "               Run step_3.sh first to compute LD scores"
            continue
        fi

        # Skip if results already exist
        if [[ -f "${OUT_PREFIX}.results" ]]; then
            log_message "      Results already exist, skipping..."
            continue
        fi

        log_message "      Running S-LDSC..."
        python "$LDSC_WRAPPER" "$LDSC_DIR" ldsc.py \
            --h2 "$SUMSTATS_FILE" \
            --ref-ld-chr "${BASELINE_LD_DIR}/baselineLD.,${CUSTOM_LD_DIR}/${REGION}." \
            --w-ld-chr "${WEIGHTS_DIR}/weights.hm3_noMHC." \
            --frqfile-chr "${FRQ_DIR}/1000G.EUR.QC." \
            --overlap-annot \
            --thin-annot \
            --print-coefficients \
            --out "$OUT_PREFIX"

        # Check if analysis succeeded
        if [[ -f "${OUT_PREFIX}.results" ]]; then
            log_message "      Success: ${OUT_PREFIX}.results"
        else
            log_message "      ERROR: S-LDSC failed for ${DISEASE}/${REGION}"
        fi
    done
done

conda deactivate
# -----------------------------------------------------------------------------
# Verify outputs
# -----------------------------------------------------------------------------
log_message "Verifying S-LDSC results..."
echo ""
echo "=== S-LDSC Results Summary ==="
echo ""

# Count results by category
neuronal_count=0
immune_count=0
vascular_count=0
control_count=0
total_count=0

for DISEASE in "${DISEASES[@]}"; do
    for REGION in "${BRAIN_REGIONS[@]}"; do
        RESULT_FILE="./results/${DISEASE}/${REGION}/${DISEASE}_${REGION}_EA.results"
        if [[ -f "$RESULT_FILE" ]]; then
            ((total_count++))
            case $DISEASE in
                ad|scz|mdd|bip|pd|smoking) ((neuronal_count++)) ;;
                ms|ra|asthma) ((immune_count++)) ;;
                cad|stroke) ((vascular_count++)) ;;
            esac
        fi
    done
done

n_regions=${#BRAIN_REGIONS[@]}
echo "Category breakdown:"
echo "  Neuronal (ad, scz, mdd, bip, pd, smoking): ${neuronal_count}/$((6 * n_regions)) results"
echo "  Immune (ms, ra, asthma): ${immune_count}/$((3 * n_regions)) results"
echo "  Vascular (cad, stroke): ${vascular_count}/$((2 * n_regions)) results"
echo ""
echo "Total: ${total_count}/${total} results"

if [[ $total_count -eq $total ]]; then
    echo ""
    echo "[OK] All S-LDSC analyses completed successfully"
else
    echo ""
    echo "[INCOMPLETE] Some analyses are missing. Check logs for errors."
fi

conda deactivate
log_message "**** Job ends ****"
