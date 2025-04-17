#!/usr/bin/env Rscript
library('methylKit')
args <- commandArgs(trailingOnly = TRUE)
my.methRaw=processBismarkAln( location = args[1],
                              sample.id=args[2], assembly="hg38", 
                              read.context="CpG", save.folder=getwd())
