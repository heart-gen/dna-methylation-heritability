#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics-gpu
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=40G

MAP_FILE="chunk_map.tsv"
> "$MAP_FILE"  # clear file

global_index=0
for CHR in {1..22}; do
    REGION_DIR="./chunked_cpg/chr_${CHR}"
    if [[ ! -d "$REGION_DIR" ]]; then
        echo "Skipping chr${CHR} (no directory)"
        continue
    fi
    for CHUNK_FILE in "$REGION_DIR"/chunk_*.bed; do
        [[ -f "$CHUNK_FILE" ]] || continue
        echo -e "${CHR}\t${CHUNK_FILE}" >> "$MAP_FILE"
        ((global_index++))
    done
done

if (( global_index == 0 )); then
    echo "No chunk files found in any chromosome directories."
    exit 1
fi

echo "Prepared mapping for $global_index jobs."

# Submit as one big array
sbatch --array=0-$((global_index - 1)) extract_snps.sh
