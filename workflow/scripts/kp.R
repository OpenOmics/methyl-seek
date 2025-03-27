#!/usr/bin/env Rscript

library(karyoploteR)
library(GenomicRanges)
library(regioneR)

args = commandArgs(trailingOnly=TRUE)

data = read.delim(args[1],h=T)
head(data)
 
D=na.omit(data)
colnames(D) = c("chrom", "start", "end", "min_p", "n_probes", "z_p", "z_sidak_p", "length")
D1 <- D[D$z_sidak_p <= 0.05, ]
D1$pval = as.numeric(as.character(D1$z_sidak_p))
D1r = toGRanges(D1)

cpgs <- createRandomRegions(nregions= D1r$length, length.mean = mean(D1r$length), mask=NA, non.overlapping = T)

for(chrom in list(c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15","chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))){
	kp <- plotKaryotype("hg19", chromosomes=chrom)
	kpPlotRegions(kp, data=cpgs, r0=0, r1=0.5)
	kpPlotDensity(kp, data=cpgs, r0=0.5, r1=1)

}
