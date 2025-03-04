
library(bsseq)
library(BiocParallel)


args=commandArgs(TRUE)

methfile <- args[1]
sid <- args[2]
cutoff <- as.numeric(args[3])
outfile <- args[4]


methylationData <- read.bismark(methfile, strandCollapse = TRUE, rmZeroCov=FALSE ,verbose=TRUE, BPPARAM =  MulticoreParam(workers = 8))
sampleNames(methylationData) <- sid
cov = getCoverage(methylationData, type="Cov", what="perBase")

sprintf("Total number of sites, %s",dim(methylationData)[1])
keepLoci.ex = which(cov >= cutoff)
meth.filtered <-  methylationData[keepLoci.ex,]
sprintf("Total number of sites after filtering low coverage (>= %s), %s",cutoff,dim(meth.filtered)[1])


meth.granges <- granges(meth.filtered)
#df.granges <- data.frame(chromosome=seqnames(meth.granges),pos=start(meth.granges))
df.granges <- data.frame(chromosome=seqnames(meth.granges),start=start(meth.granges),end=end(meth.granges)+1)



df <- data.frame(df.granges,data.frame(getMeth(meth.filtered, type="raw")))
write.csv(df,outfile,quote = FALSE,row.names = FALSE)


print(sprintf("%s, done",sid))


