#### Perform principal component analysis ####

library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')
library(scales)

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