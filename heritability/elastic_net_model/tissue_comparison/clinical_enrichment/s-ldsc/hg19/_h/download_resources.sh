#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --job-name=download_ldsc_resources
#SBATCH --output=logs/download_resources.output_%j.log
#SBATCH --error=logs/download_resources.error_%j.log

# =============================================================================
# Download LDSC Reference Files from Broad Institute
# =============================================================================
# This script downloads the necessary reference files for running S-LDSC analysis
# Files are downloaded to /ocean/projects/bio250020p/shared/resources/ldsc/
# =============================================================================

set -e

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

echo "**** PSC Bridges-2 info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"

# Resource directory
RESOURCE_DIR="/ocean/projects/bio250020p/shared/resources/ldsc"
mkdir -p "$RESOURCE_DIR"
cd "$RESOURCE_DIR"

# URLs for LDSC resources from Broad Institute
BASE_URL="https://data.broadinstitute.org/alkesgroup/LDSCORE"

# Resources to download
declare -A RESOURCES=(
    ["1000G_Phase3_plinkfiles"]="${BASE_URL}/1000G_Phase3_plinkfiles.tgz"
    ["1000G_Phase3_baselineLD_v2.2_ldscores"]="${BASE_URL}/1000G_Phase3_baselineLD_v2.2_ldscores.tgz"
    ["1000G_Phase3_weights_hm3_no_MHC"]="${BASE_URL}/1000G_Phase3_weights_hm3_no_MHC.tgz"
    ["1000G_Phase3_frq"]="${BASE_URL}/1000G_Phase3_frq.tgz"
)

# Download and extract each resource
for RESOURCE_NAME in "${!RESOURCES[@]}"; do
    URL="${RESOURCES[$RESOURCE_NAME]}"
    FILENAME=$(basename "$URL")

    log_message "Processing: $RESOURCE_NAME"

    # Check if already extracted
    if [[ -d "$RESOURCE_NAME" ]]; then
        log_message "  Directory $RESOURCE_NAME already exists, skipping..."
        continue
    fi

    # Download if not present
    if [[ ! -f "$FILENAME" ]]; then
        log_message "  Downloading $FILENAME..."
        wget -q --show-progress "$URL" -O "$FILENAME"
    else
        log_message "  Archive $FILENAME already exists, skipping download..."
    fi

    # Extract
    log_message "  Extracting $FILENAME..."
    tar -xzf "$FILENAME"

    # Clean up archive to save space (optional - comment out to keep archives)
    # rm "$FILENAME"

    log_message "  Done with $RESOURCE_NAME"
done

# Download HapMap3 SNP list if not present
HAPMAP3_FILE="w_hm3.snplist"
HAPMAP3_URL="${BASE_URL}/${HAPMAP3_FILE}.bz2"

if [[ ! -f "$HAPMAP3_FILE" ]]; then
    log_message "Downloading HapMap3 SNP list..."
    wget -q --show-progress "$HAPMAP3_URL" -O "${HAPMAP3_FILE}.bz2"
    bunzip2 "${HAPMAP3_FILE}.bz2"
fi

# Download liftover chain file
LIFTOVER_DIR="/ocean/projects/bio250020p/shared/resources/liftover"
mkdir -p "$LIFTOVER_DIR"

CHAIN_FILE="${LIFTOVER_DIR}/hg38ToHg19.over.chain.gz"
CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"

if [[ ! -f "$CHAIN_FILE" ]]; then
    log_message "Downloading hg38ToHg19 liftover chain file..."
    wget -q --show-progress "$CHAIN_URL" -O "$CHAIN_FILE"
else
    log_message "Chain file already exists: $CHAIN_FILE"
fi

# Verify downloads
log_message "Verifying downloads..."
echo ""
echo "=== Downloaded Resources ==="
for RESOURCE_NAME in "${!RESOURCES[@]}"; do
    if [[ -d "$RESOURCE_DIR/$RESOURCE_NAME" ]]; then
        echo "[OK] $RESOURCE_NAME"
    else
        echo "[MISSING] $RESOURCE_NAME"
    fi
done

if [[ -f "$CHAIN_FILE" ]]; then
    echo "[OK] Liftover chain file"
else
    echo "[MISSING] Liftover chain file"
fi

echo ""
log_message "**** Job ends ****"
