#### Extract variably methylated regions ####

library('bsseq')
library('HDF5Array')
library(DelayedMatrixStats)
library('data.table')
library(scales)
library(GenomicRanges)

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

#extract methylation values for 1 vmr
chr <- "chr1"
start <- 903969
end <- 904084
reg <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
meth_reg <- getMeth(BSobj, regions = reg, what="perRegion")