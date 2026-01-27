library(GenomicRanges)
library(rtracklayer)

library(rtracklayer)

# Load BED file
gr <- import("gwas.bed", format="bed")

# Load chain
chain <- import.chain("hg19ToHg38.over.chain")

# LiftOver
lifted <- liftOver(gr, chain)
lifted_gr <- unlist(lifted)

# Export
export(lifted_gr, "gwas_lifted.bed", format="bed")

