#### VMR Test ####

library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')
library(scales)

### Adapting sd.R 01_pca.R and 03_vmr.R for one chromosome
load_bsobj <- function(working_dir, file_name) {
    setwd(working_dir)
    load(file_name)
}

# Change file path for raw data
update_bsobj <- function(BSobj, raw_path) {
    var <- c("M","Cov","coef")
    for (i in var) {
        BSobj@assays@data@listData[[i]]@seed@seed@filepath <- raw_path
    }
    return(BSobj)
}
#BSobj = BSobj[seqnames(BSobj) == "chr1",] #keep chr1

# keep only control, adult AA
pheno <- "/projects/p32505/projects/dna-methylation-heritability/inputs/phenotypes/merged/_m/merged_phenotypes.csv" 
ances <- read.csv(pheno,header=T)
id <- intersect(ances$BrNum[ances$Race == "AA" & ances$Dx == "Control" & ances$Age >= 17],colData(BSobj)$brnum)
BSobj <- BSobj[,is.element(colData(BSobj)$brnum,id)]

# exlcude low coverage sites
cov=getCoverage(BSobj)
n <- length(colData(BSobj)$brnum)
keep <- which(rowSums2(cov >= 5) >= n * 0.8) 
BSobj <- BSobj[keep,]

# calculate sd & mean of DNAm
M=as.matrix(getMeth(BSobj))
sds <- rowSds(M)
means <- rowMeans2(M)
save(sds,means,BSobj,file=paste0("./","chr1_stats",".rda"))

# load sd of raw DNAm
i <- 1
v <- data.frame(chr=i,start=start(BSobj),sd=sds)

# get top 1M variable CpG
v <- v[order(v$sd * -1),]
v <- v[1:10^6,]

# DNAm levels for top 1M
tmp <- v[v$chr==i,]
BS <- BSobj[is.element(start(BSobj),tmp$start),]
meth <- as.matrix(getMeth(BS))

id_top <- colData(BS)$brnum
colnames(meth) <- id_top

# pca
pc = prcomp(t(meth), scale. = T, center = T)
cumsum(summary(pc)$importance[2,1:10])
# PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10
#0.11994 0.16300 0.19064 0.21188 0.22916 0.24498 0.25889 0.27264 0.28560 0.29838

pdf("pca.pdf")
plot(pc)
plot(pc$x[,1], pc$x[,2])
dev.off()
write.csv(pc$x,"pc.csv")

# read files
#fs <- list.files(path = ".", pattern = "^p",recursive = T, full.names = T)
#fs <- fs[!grepl("pc",fs)]
#raw <-  lapply(fs, function(x) fread(x,data.table=F))
#v <- rbindlist(raw)
#colnames(v) <- c("chr","start","sd")

# sd cutoff
sdCut <- quantile(v[,3], prob = 0.99, na.rm = TRUE)

# get vmr
vmrs <- c()
#	v2 <- v[v$chr==paste0("chr",i),] # Error
#	v2 <- v2[order(v2$start)]
isHigh <- rep(0, nrow(v))
isHigh[v$sd > sdCut] <- 1
vmrs0 <- bsseq:::regionFinder3(isHigh, as.character(v$chr), v$start, maxGap = 1000)$up
vmr <- vmrs0[vmrs0$n > 5,1:3]
write.table(vmr,"vmr.bed",col.names=F,row.names=F,sep="\t",quote=F)

# adapt for all chr
#for(i in 1:1){
#	cat(i,"\n")
#	v2 <- v[v$chr==paste0("chr",i),]
#	v2 <- v2[order(v2$start)]
#	isHigh <- rep(0, nrow(v2))
#	isHigh[v2$sd > sdCut] <- 1
#	vmrs0 <- bsseq:::regionFinder3(isHigh, as.character(v2$chr), v2$start, maxGap = 1000)$up
#	vmrs <- rbind(vmrs,vmrs0)
#}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()