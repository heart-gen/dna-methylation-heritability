#### Stratify SNPs by LD score of individual variants ####

suppressPackageStartupMessages({
    library(here)
    library(ggplot2)
})

args  <- commandArgs(trailingOnly = TRUE)
chr   <- args[1]

## Function
get_snp_groups <- function(lds_seg) {
    quartiles <- summary(lds_seg$ldscore_SNP)
    snp_groups <- list(
        lb1 = which(lds_seg$ldscore_SNP <= quartiles[2]),
        lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3]),
        lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5]),
        lb4 = which(lds_seg$ldscore_SNP > quartiles[5])
    )
    return(list(quartiles = quartiles, snp_groups = snp_groups))
}

extract_snps <- function(lds_seg, snp_groups, out_ld) {
    snps <- list(
      lb1_snp = lds_seg$SNP[snp_groups$lb1],
      lb2_snp = lds_seg$SNP[snp_groups$lb2],
      lb3_snp = lds_seg$SNP[snp_groups$lb3],
      lb4_snp = lds_seg$SNP[snp_groups$lb4]
    )
    
    write.table(snps$lb1_snp, file=file.path(out_ld,"snp_group_1.txt"), row.names=F, quote=F, col.names=F)
    write.table(snps$lb2_snp, file=file.path(out_ld,"snp_group_2.txt"), row.names=F, quote=F, col.names=F)
    write.table(snps$lb3_snp, file=file.path(out_ld,"snp_group_3.txt"), row.names=F, quote=F, col.names=F)
    write.table(snps$lb4_snp, file=file.path(out_ld,"snp_group_4.txt"), row.names=F, quote=F, col.names=F)
    return(snps)
}

plot_ld <- function(lds_seg, quartiles, chr, out_ld) {
    hist <- ggplot(lds_seg, aes(x = ldscore_SNP)) + 
        geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black", alpha = 0.7) +
        geom_vline(aes(xintercept = quartiles[2]), color = "red", linetype = "dashed") +
        geom_vline(aes(xintercept = quartiles[3]), color = "red", linetype = "dashed") +  
        geom_vline(aes(xintercept = quartiles[5]), color = "red", linetype = "dashed") +  
        labs(title = paste("Chr", chr, "SNP LD Scores"),
            x = "SNP LD Score", y = "Frequency")

    pdf(file=file.path(out_ld,"ld_hist.pdf"))
    print(hist)
    dev.off()
    return(hist)
}

## Main
# Load LD scores generated from gcta
lds_seg = read.table(
  here("testing/caudate/_m/h2", 
      paste0("chr_", chr, "/TOPMed_LIBD.AA.VMR1.score.ld")),
  header = TRUE,
  colClasses = c("character", rep("numeric", 8))
)
out_ld <- here("testing/caudate/_m/h2", paste0("chr_", chr))

# Stratify snps based on LD scores
groups <- get_snp_groups(lds_seg)

# Extract each group of snps and write to text file
snps <- extract_snps(lds_seg, groups$snp_groups, out_ld)

# Plot histogram 
hist <- plot_ld(lds_seg, groups$quartiles, chr, out_ld)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()