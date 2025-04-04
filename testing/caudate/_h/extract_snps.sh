#!/bin/bash

#module load plink/2.0-alpha-3.3

data_dir="/projects/p32505/projects/dna-methylation-heritability/inputs/genotypes"
output_dir="/projects/p32505/users/alexis/projects/dna-methylation-heritability/testing/caudate/_m/chr_1"
chr_file="/projects/b1213/resources/genomes/human/gencode-v47/fasta/chromosome_sizes.txt"
ranges="/projects/p32505/users/alexis/projects/dna-methylation-heritability/testing/caudate/_m/chr_1/vmr.bed"

#for chr in {1..22}; do
#mkdir /projects/p32505/users/alexis/projects/dna-methylation-heritability/testing/caudate/_m/chr_$chr ;
#done

while IFS=$'\t' read -r chrom start end; do 
    echo "Processing Chromosome $chrom: $start-$end"

    plink2 --pfile $data_dir/TOPMed_LIBD.AA --chr "$chrom" --from-bp "$start" --to-bp "$end" --make-bed --out "$output_dir/TOPMed_LIBD.AA.Chr${chrom}_${start}_${end}"

    # check chromosome size information
    window=500000
    chrom_size=$(grep "^chr1[[:space:]]" $chr_file | cut -f2)

    start_pos=$start-$window

    if (( start_pos < 1 )); then
        echo "ERROR: Start position is below zero."
    fi

    end_pos=$end+$window

    if (( end_pos >= chrom_size )); then
        echo "ERROR: End position exceeds Chromosome $chrom size."
    fi

done < "$ranges"

# extract snps for 1 vmr
# plink2 --pfile $data_dir/TOPMed_LIBD.AA --from-bp 403969 --to-bp 1404084 --maf 0.05 --geno 0.1 --make-bed --out $output_dir/TOPMed_LIBD.AA.VMR1


