#!/bin/bash

#module load gcta/1.94.0
#module load R/4.2
##### GREML-LDMS #####
gcta64 --bfile TOPMed_LIBD.AA.VMR1 --ld-score-region 200 --out TOPMed_LIBD.AA.VMR1

Rscript 05.stratify_LD.R

# Make GRM for each group
gcta64 --bfile TOPMed_LIBD.AA.VMR1 --extract snp_group1.txt --make-grm --out TOPMed_LIBD.AA.VMR1_group1
gcta64 --bfile TOPMed_LIBD.AA.VMR1 --extract snp_group2.txt --make-grm --out TOPMed_LIBD.AA.VMR1_group2
gcta64 --bfile TOPMed_LIBD.AA.VMR1 --extract snp_group3.txt --make-grm --out TOPMed_LIBD.AA.VMR1_group3
gcta64 --bfile TOPMed_LIBD.AA.VMR1 --extract snp_group4.txt --make-grm --out TOPMed_LIBD.AA.VMR1_group4

# GREML with multiple GRM
gcta64 --reml --mgrm multi_GRMs.txt --pheno methylation.phen --out TOPMed_LIBD.AA.VMR1 #Error: analysis stopped because more than half of the variance components are constrained. The result would be unreliable.

# REML no contrain 
gcta64 --reml-no-constrain --mgrm multi_GRMs.txt --pheno methylation.phen --out TOPMed_LIBD.AA.VMR1 #Error: the information matrix is not invertible.
