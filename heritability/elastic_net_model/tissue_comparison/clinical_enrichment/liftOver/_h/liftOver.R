library(GenomicRanges)
library(rtracklayer)

library(rtracklayer)

# Load BED file
gr <- import("AFR.bed", format="bed")

# Load chain
chain <- import.chain("hg19ToHg38.over.chain")

# LiftOver
lifted <- liftOver(gr, chain)
lifted_gr <- unlist(lifted)

# Export
export(lifted_gr, "AFR_lifted.bed", format="bed")

