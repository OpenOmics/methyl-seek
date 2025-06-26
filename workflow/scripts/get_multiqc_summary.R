#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(stringr)

aln_file = args[1]
dup_file = args[2]
meth_file = args[3]
fas_file = args[4]

aln = read.delim(aln_file)
aln$Library_ID = str_split_fixed(aln$Sample, '[.]', 2)[,1]
aln1 = aln[,c("Library_ID","aligned_reads","percent_aligned")]

dup = read.delim(dup_file)
dup$Library_ID = str_split_fixed(dup$Sample, '[.]', 3)[,1]
dup1 = dup[,c("Library_ID","dup_reads_percent")]

meth = read.delim(meth_file)
meth$Library_ID = str_split_fixed(meth$Sample, '[.]', 2)[,1]
meth1 = meth[,c("Library_ID","meth_cpg")]

fas = read.delim(meth_file)
fas$Library_ID = str_split_fixed(fas$Sample, '[.]', 2)[,1]

fas1 = aggregate(fas$Total.Sequences, list(fas$Library_ID), FUN=mean)
fas2 = aggregate(fas$Sequences.flagged.as.poor.quality, list(fas$Library_ID), FUN=mean)
fas1$y = fas1$x-fas2$x
colnames(fas1) = c("Library_ID","Total.pairs","Filtered.pairs")

total = merge(merge(merge(x=fas1,y=aln1,by="Library_ID",all=T),dup1,by="Library_ID",all=T),meth1,by="Library_ID",all=T)
total$BS_conv_rate = total$meth_cpg/total$aligned_reads

total = total[,-7]

write.csv(x = total,file = args[5],row.names = F,quote = F)
