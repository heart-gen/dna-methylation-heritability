#!/bin/bash

# === CONFIG ===
total_tasks=9000
chunk_size=750
script="../_h/elastic_h2_array-01.sh"
jobids=()

# === SUBMIT JOBS IN CHUNKS WITH DEPENDENCY ===
for ((offset=4501; offset<=total_tasks; offset+=chunk_size)); do
    end=$((offset + chunk_size - 1))
    [ "$end" -gt "$total_tasks" ] && end="$total_tasks"
    count=$((end - offset + 1))

    echo "Submitting array for OFFSET=$offset (tasks $offset to $end)"

    jobid=$(sbatch --account=p32505 --partition=short --mem=10GB --array=1-${count}%250 --time=00:30:00 --ntasks=1 --cpus-per-task=1 --job-name=elastic_h2 --mail-type=FAIL --mail-user=alexis.bennett@northwestern.edu --output=logs/elastic_h2_%A_%a.log --export=OFFSET=${offset},ALL --parsable "$script")
    echo "Submitted JobID: $jobid for OFFSET=$offset"
    jobids+=("$jobid")
done
