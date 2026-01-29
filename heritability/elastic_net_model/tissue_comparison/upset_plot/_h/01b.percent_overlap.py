"""
This script calculates the percent overlap of VMRs
across brain regions
"""

import session_info

input_fn = "./F_0.5/sets/caudate_hippocampus_overlap.bed" 
output_fn = "./F_0.5/percent_overlap/caudate_hippocampus_overlap.tsv" 

with open(input_fn) as f, open(output_fn, "w") as out:
    out.write("chrom\tstartA\tendA\tstartB\tendB\toverlap\tpctA\tpctB\treciprocal_overlap\n")
    for line in f:
        cols = line.strip().split()
        
        start_a, end_a = int(cols[1]), int(cols[2])
        start_b, end_b = int(cols[4]), int(cols[5])
        overlap = int(cols[6])

        len_a = end_a - start_a
        len_b = end_b - start_b

        # Percent overlap for each VMR
        pct_a = (overlap / len_a) * 100
        pct_b = (overlap / len_b) * 100

        # Reciprocal overlap
        reciprocal = min(pct_a, pct_b)

        out.write(f"{cols[0]}\t{start_a}\t{end_a}\t{start_b}\t{end_b}\t{overlap}\t"
                  f"{pct_a:.2f}\t{pct_b:.2f}\t{reciprocal:.2f}\n")

# Session information
session_info.show()