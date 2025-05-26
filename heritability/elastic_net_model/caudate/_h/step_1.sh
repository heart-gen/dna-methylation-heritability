#!/bin/bash

total_jobs=12078
max_chunk=1000
script="../_h/elastic_h2_array-01.sh"

mkdir -p logs

for ((start=1; start<=total_jobs; start+=max_chunk)); do
    end=$((start + max_chunk - 1))
    if [ $end -gt $total_jobs ]; then
        end=$total_jobs
    fi
    echo "Submitting array job: $start-$end"
    sbatch --array=$start-$end $script
done
