#### Perform principal component analysis ####

library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')

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

# get ancestry
f_ances <- here("inputs/genetic-ancestry/structure.out_ancestry_proportion_raceDemo_compare")
gen_ances <- read.table(f_ances,header=T)
idx <- match(id_top,gen_ances$id)

# residual after regressing out Afr proportion
res <- meth
for(i in 1:nrow(meth)){
  if(! i %% 10000){
    cat(i,"\n")
  }
  model = lm(meth[i,] ~ gen_ances$Afr[idx])
  res[i,] = resid(model)
}

# pca
pc = prcomp(t(res), scale. = T, center = T)
cumsum(summary(pc)$importance[2,1:10])
#PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10 
#0.05648 0.08545 0.11074 0.13106 0.15105 0.16687 0.17921 0.19071 0.20140 0.21163 

pdf(file=file.path(output,"pca.pdf"))
plot(pc)
plot(pc$x[,1], pc$x[,2])
dev.off()
write.csv(pc$x,file=file.path(output,"pc.csv"))

# correlation of pc with ancestry
res <- matrix(NA,nrow=ncol(pc$x),ncol=2)
for(i in 1:ncol(pc$x)){
  tmp <- cor.test(gen_ances$Afr[idx],pc$x[,i])
  res[i,1] <- tmp$estimate
  res[i,2] <- tmp$p.value
}
write.csv(res,file=file.path(output,"pc_ances_cor.csv"))