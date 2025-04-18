#!/bin/bash

#module load gcta/1.94.0
#module load R/4.2
##### GREML-LDMS #####
gcta64 --bfile TOPMed_LIBD.AA.VMR1 --ld-score-region 200 --out TOPMed_LIBD.AA.VMR1

Rscript 05.stratify_LD.R

# Make GRM for each group
for i in 1:4 ; do
gcta64 --bfile TOPMed_LIBD.AA.VMR1 --extract snp_group_${i}.txt --make-grm --out TOPMed_LIBD.AA.VMR1_group_${i}

echo "TOPMed_LIBD.AA.VMR1_group_${i}" >> VMR1_multi_GRMs.txt
done

# GREML with multiple GRM
gcta64 --reml --mgrm VMR1_multi_GRMs.txt --pheno VMR1_meth.phen --covar TOPMed_LIBD.AA.covar --qcovar TOPMed_LIBD.AA.qcovar --out TOPMed_LIBD.AA.VMR1