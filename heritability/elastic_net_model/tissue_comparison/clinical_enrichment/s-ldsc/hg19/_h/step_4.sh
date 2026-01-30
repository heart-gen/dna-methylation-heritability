#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --job-name=partition_heritability
#SBATCH --output=logs/04.output_%j.log
#SBATCH --error=logs/04.error_%j.log

# =============================================================================
# Step 4: Partition Heritability with S-LDSC
# =============================================================================
# This script runs stratified LD score regression (S-LDSC) to partition
# heritability for each disease/trait across brain regions and heritability
# status. Processes 10 diseases x 3 regions x 2 heritability = 60 analyses.
# =============================================================================

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

# Validate resources
if ! validate_resources; then
    log_message "ERROR: Resource validation failed. Exiting."
    exit 1
fi

# LDSC script
LDSC_SCRIPT="${LDSC_DIR}/ldsc.py"

# Resource paths
BASELINE_LD_DIR="${RESOURCE_DIR}/1000G_Phase3_baselineLD_v2.2_ldscores"
WEIGHTS_DIR="${RESOURCE_DIR}/1000G_Phase3_weights_hm3_no_MHC"
FRQ_DIR="${RESOURCE_DIR}/1000G_Phase3_frq"

# Summary statistics directory
SUMSTATS_DIR="./sumstats"

# 10 diseases spanning neuronal, immune, and vascular categories
DISEASES=("ad" "scz" "mdd" "bip" "pd" "ms" "ra" "asthma" "cad" "htn")

# Counter for progress
total=$((${#DISEASES[@]} * ${#BRAIN_REGIONS[@]} * ${#HERITABILITY[@]}))
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
        log_message "  Processing region: ${REGION}"

        # Process each heritability status
        for STATUS in "${HERITABILITY[@]}"; do
            ((current++))
            log_message "    Processing status: ${STATUS} (${current}/${total})"

            # Custom LD scores directory
            CUSTOM_LD_DIR="./custom_ldscores/${REGION}/${STATUS}"

            # Output directory
            OUT_DIR="./results/${DISEASE}/${REGION}/${STATUS}"
            mkdir -p "$OUT_DIR"

            # Output file prefix
            OUT_PREFIX="${OUT_DIR}/${DISEASE}_${REGION}_${STATUS}"

            # Check if custom LD scores exist
            if [[ ! -f "${CUSTOM_LD_DIR}/${REGION}_${STATUS}.1.l2.ldscore.gz" ]]; then
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
            python "$LDSC_SCRIPT" \
                --h2 "$SUMSTATS_FILE" \
                --ref-ld-chr "${BASELINE_LD_DIR}/baselineLD.,${CUSTOM_LD_DIR}/${REGION}_${STATUS}." \
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
                log_message "      ERROR: S-LDSC failed for ${DISEASE}/${REGION}/${STATUS}"
            fi
        done
    done
done

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
total_count=0

for DISEASE in "${DISEASES[@]}"; do
    for REGION in "${BRAIN_REGIONS[@]}"; do
        for STATUS in "${HERITABILITY[@]}"; do
            RESULT_FILE="./results/${DISEASE}/${REGION}/${STATUS}/${DISEASE}_${REGION}_${STATUS}.results"
            if [[ -f "$RESULT_FILE" ]]; then
                ((total_count++))
                case $DISEASE in
                    ad|scz|mdd|bip|pd) ((neuronal_count++)) ;;
                    ms|ra|asthma) ((immune_count++)) ;;
                    cad|htn) ((vascular_count++)) ;;
                esac
            fi
        done
    done
done

echo "Category breakdown:"
echo "  Neuronal (ad, scz, mdd, bip, pd): ${neuronal_count}/30 results"
echo "  Immune (ms, ra, asthma): ${immune_count}/18 results"
echo "  Vascular (cad, htn): ${vascular_count}/12 results"
echo ""
echo "Total: ${total_count}/60 results"

if [[ $total_count -eq 60 ]]; then
    echo ""
    echo "[OK] All S-LDSC analyses completed successfully"
else
    echo ""
    echo "[INCOMPLETE] Some analyses are missing. Check logs for errors."
fi

log_message "**** Job ends ****"
