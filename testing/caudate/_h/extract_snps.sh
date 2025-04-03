#!/bin/bash

#module load plink/2.0-alpha-3.3

data_dir="/projects/p32505/projects/dna-methylation-heritability/inputs/genotypes"
output_dir="/projects/p32505/users/alexis/projects/dna-methylation-heritability/testing/caudate/_m"
chr_file="/projects/b1213/resources/genomes/human/gencode-v47/fasta/chromosome_sizes.txt"
ranges="/projects/p32505/users/alexis/projects/dna-methylation-heritability/testing/caudate/_h/ranges.txt"

# extract snps for 1 vmr
plink2 --pfile $data_dir/TOPMed_LIBD.AA --from-bp 403969 --to-bp 1404084 --maf 0.05 --geno 0.1 --make-bed --out $output_dir/TOPMed_LIBD.AA.VMR1

# check chromosome size information
window=500000
chrom=1
start_pos=403969
end_pos=1404084

chrom_size=$(grep "^chr1[[:space:]]" $chr_file | cut -f2)

start=$start_pos-$window

if (( start < 1 )); then
    echo "ERROR: Start position is below zero."
fi

end=$end_pos+$window

if (( end >= chrom_size )); then
    echo "ERROR: End position exceeds Chromosome $chrom size."
fi
