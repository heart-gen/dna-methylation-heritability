#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=make_plots
#SBATCH --output=logs/05.output_make_plots.log
#SBATCH --error=logs/05.error_make_plots.log
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# =============================================================================
# Step 5: Combine Results and Generate Plots
# =============================================================================
# This script collects all S-LDSC results, combines them into a single table,
# and generates visualization plots (heatmaps) for comparing enrichment across
# diseases and brain regions.
# =============================================================================

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

log_message "**** Job starts ****"
print_job_info

# 10 diseases spanning neuronal, immune, and vascular categories
DISEASES=("ad" "scz" "mdd" "bip" "pd" "ms" "ra" "asthma" "cad" "htn")

# Brain regions and heritability status
BRAIN_REGIONS=("caudate" "dlpfc" "hippocampus")
HERITABILITY=("heritable_hg19" "non_heritable_hg19")

# -----------------------------------------------------------------------------
# Dynamically generate file list
# -----------------------------------------------------------------------------
log_message "Generating file list..."

LDSC_FILES=()
missing_count=0
found_count=0

for DISEASE in "${DISEASES[@]}"; do
    for REGION in "${BRAIN_REGIONS[@]}"; do
        for STATUS in "${HERITABILITY[@]}"; do
            FILE="results/${DISEASE}/${REGION}/${STATUS}/${DISEASE}_${REGION}_${STATUS}.results"
            if [[ -f "$FILE" ]]; then
                LDSC_FILES+=("$FILE")
                ((found_count++))
            else
                log_message "  Missing: $FILE"
                ((missing_count++))
            fi
        done
    done
done

log_message "Found ${found_count} result files, missing ${missing_count}"

# Check if we have any files to process
if [[ ${#LDSC_FILES[@]} -eq 0 ]]; then
    log_message "ERROR: No result files found. Run step_4.sh first."
    exit 1
fi

# Write file list
printf '%s\n' "${LDSC_FILES[@]}" > ldsc_file_list.txt
log_message "File list written to ldsc_file_list.txt"

# -----------------------------------------------------------------------------
# Run Python plotting script
# -----------------------------------------------------------------------------
log_message "Running plotting script..."

# Check if Python script exists
PLOT_SCRIPT="${SCRIPT_DIR}/05.make_plots.py"
if [[ ! -f "$PLOT_SCRIPT" ]]; then
    log_message "WARNING: Plotting script not found: $PLOT_SCRIPT"
    log_message "Creating combined results table only..."

    # Create combined results table manually
    OUTPUT_FILE="combined_ldsc_results.tsv"

    # Header
    echo -e "Disease\tDisease_Name\tCategory\tRegion\tHeritability_Status\tCategory_Name\tProp_SNPs\tProp_h2\tProp_h2_std_error\tEnrichment\tEnrichment_std_error\tEnrichment_p" > "$OUTPUT_FILE"

    # Disease name mapping
    declare -A DISEASE_NAMES=(
        ["ad"]="Alzheimer's Disease"
        ["scz"]="Schizophrenia"
        ["mdd"]="Major Depression"
        ["bip"]="Bipolar Disorder"
        ["pd"]="Parkinson's Disease"
        ["ms"]="Multiple Sclerosis"
        ["ra"]="Rheumatoid Arthritis"
        ["asthma"]="Asthma"
        ["cad"]="Coronary Artery Disease"
        ["htn"]="Hypertension"
    )

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

    # Process each result file
    for FILE in "${LDSC_FILES[@]}"; do
        # Parse file path to get metadata
        # Format: results/DISEASE/REGION/STATUS/DISEASE_REGION_STATUS.results
        DISEASE=$(echo "$FILE" | cut -d'/' -f2)
        REGION=$(echo "$FILE" | cut -d'/' -f3)
        STATUS=$(echo "$FILE" | cut -d'/' -f4)

        DISEASE_NAME="${DISEASE_NAMES[$DISEASE]}"
        CATEGORY="${DISEASE_CATEGORIES[$DISEASE]}"

        # Extract the annotation row (last row, which is our custom annotation)
        # Skip header and baseline annotations
        ANNOT_ROW=$(tail -1 "$FILE")

        # Parse columns (assuming standard LDSC output format)
        # Columns: Category, Prop._SNPs, Prop._h2, Prop._h2_std_error, Enrichment, Enrichment_std_error, Enrichment_p
        CATEGORY_NAME=$(echo "$ANNOT_ROW" | awk '{print $1}')
        PROP_SNPS=$(echo "$ANNOT_ROW" | awk '{print $2}')
        PROP_H2=$(echo "$ANNOT_ROW" | awk '{print $3}')
        PROP_H2_SE=$(echo "$ANNOT_ROW" | awk '{print $4}')
        ENRICHMENT=$(echo "$ANNOT_ROW" | awk '{print $5}')
        ENRICHMENT_SE=$(echo "$ANNOT_ROW" | awk '{print $6}')
        ENRICHMENT_P=$(echo "$ANNOT_ROW" | awk '{print $7}')

        echo -e "${DISEASE}\t${DISEASE_NAME}\t${CATEGORY}\t${REGION}\t${STATUS}\t${CATEGORY_NAME}\t${PROP_SNPS}\t${PROP_H2}\t${PROP_H2_SE}\t${ENRICHMENT}\t${ENRICHMENT_SE}\t${ENRICHMENT_P}" >> "$OUTPUT_FILE"
    done

    log_message "Combined results written to: $OUTPUT_FILE"
else
    python "$PLOT_SCRIPT"
fi

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
log_message "Generating summary..."
echo ""
echo "=== S-LDSC Analysis Summary ==="
echo ""
echo "Diseases analyzed: ${#DISEASES[@]}"
echo "  Neuronal: ad, scz, mdd, bip, pd"
echo "  Immune: ms, ra, asthma"
echo "  Vascular: cad, htn"
echo ""
echo "Brain regions: ${BRAIN_REGIONS[*]}"
echo "Heritability status: ${HERITABILITY[*]}"
echo ""
echo "Total analyses: $((${#DISEASES[@]} * ${#BRAIN_REGIONS[@]} * ${#HERITABILITY[@]}))"
echo "Completed: ${found_count}"
echo "Missing: ${missing_count}"
echo ""

if [[ -f "combined_ldsc_results.tsv" ]]; then
    echo "Output files:"
    echo "  - combined_ldsc_results.tsv"
    echo "  - ldsc_file_list.txt"
    if [[ -f "enrichment_heatmap.pdf" ]]; then
        echo "  - enrichment_heatmap.pdf"
    fi
fi

log_message "**** Job ends ****"
