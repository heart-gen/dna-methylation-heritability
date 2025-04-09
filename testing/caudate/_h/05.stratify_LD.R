setwd("/projects/p32505/projects/dna-methylation-heritability/testing/_h")
lds_seg = read.table("../_m/TOPMed_LIBD.AA.VMR1.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)

lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]

hist <- ggplot(lds_seg, aes(x = ldscore_SNP)) + 
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = quartiles[2]), color = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = quartiles[3]), color = "red", linetype = "dashed") +  
  geom_vline(aes(xintercept = quartiles[5]), color = "red", linetype = "dashed") +  
  labs(title = "Chr 1 SNP LD Scores",
       x = "SNP LD Score", y = "Frequency")

pdf(file=file.path(output,"ld_hist.pdf"))
print(hist)
dev.off()

write.table(lb1_snp, file=file.path(output,"snp_group_1.txt"), row.names=F, quote=F, col.names=F)
write.table(lb2_snp, file=file.path(output,"snp_group_2.txt"), row.names=F, quote=F, col.names=F)
write.table(lb3_snp, file=file.path(output,"snp_group_3.txt"), row.names=F, quote=F, col.names=F)
write.table(lb4_snp, file=file.path(output,"snp_group_4.txt"), row.names=F, quote=F, col.names=F)
