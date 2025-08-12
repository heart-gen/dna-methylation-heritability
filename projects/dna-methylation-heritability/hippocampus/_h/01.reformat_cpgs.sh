#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1G
#SBATCH --job-name=reformat_cpgs
#SBATCH --output=logs/reformat_cpgs_out.log
#SBATCH --error=logs/reformat_cpgs_err.log

mkdir -p logs

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

echo "**** QUEST info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

input_file_name="cpg_pos.txt"
chrom_sizes_file="/projects/b1213/resources/genomes/human/gencode-v47/fasta/chromosome_sizes.txt"

if [[ -z "$input_file_name" ]] || [[ ! -f "$chrom_sizes_file" ]]; then
    log_message "Usage: $0 — missing $input_file_name or $chrom_sizes_file"
    exit 1
fi

# Load chromosome sizes, stripping "chr" prefix for array keys
declare -A chrom_sizes
while read -r chrom size; do
    chrom_no="${chrom#chr}"   # remove 'chr' prefix
    chrom_sizes["$chrom_no"]=$size
done < "$chrom_sizes_file"

cd cpg

# Find dirs like chr_1, chr_X, chr_Y ...
find . -type d -regex ".*/chr_[0-9XY]+" | while read -r dir; do
    input_file="$dir/$input_file_name"

    if [[ ! -f "$input_file" ]]; then
        log_message "[SKIP] $dir — $input_file_name not found."
        continue
    fi

    chr_name=$(basename "$dir")
    if [[ "$chr_name" =~ ^chr_([0-9XY]+)$ ]]; then
        chrom="${BASH_REMATCH[1]}"  # no "chr" prefix
    else
        log_message "[WARN] Invalid directory name: $dir"
        continue
    fi

    chrom_size=${chrom_sizes["$chrom"]}
    if [[ -z "$chrom_size" ]]; then
        log_message "[WARN] Chromosome size not found for $chrom — skipping"
        continue
    fi

    log_message "[INFO] Processing chromosome $chrom ($input_file)..."

    output_file="$dir/formatted_cpg.bed"
    row=1
    site_count=0
    skipped_count=0

    > "$output_file"

    while IFS= read -r X; do
        bed_start=$((X - 20000))      # no -1, keep 1-based
        bed_end=$((X + 20000))        # inclusive end

        # Skip if out of bounds
        if (( bed_start <= 0 || bed_end >= chrom_size )); then
            skipped_count=$((skipped_count + 1))
            continue
        fi

        echo -e "${chrom}\t${bed_start}\t${bed_end}\tCpG_${row}" >> "$output_file"
        row=$((row + 1))
        site_count=$((site_count + 1))
    done < <(tail -n +3 "$input_file")

    log_message "[DONE] Wrote $site_count CpG sites to $output_file (skipped $skipped_count out-of-bounds)"
done

log_message "**** Job ends ****"
